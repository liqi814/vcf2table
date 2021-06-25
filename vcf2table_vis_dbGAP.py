import pandas as pd
import io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


with open('dbGAP/data/CaseSampleName_ID_dict_hg38.txt') as dict_file:
    CaseQV_pair = dict(x.rstrip().split(None, 1) for x in dict_file)
phased_VCF = read_vcf('dbGAP/data/dbGAP_IPF.KIF15.ID_nomono.vcf')
with open('dbGAP/data/KIF15_topsnps_hg38_shortlist.txt', 'r') as f:
    common_snps = f.readlines()
# remove whitespace characters like `\n` at the end of each line
common_snps = [x.strip() for x in common_snps]

phased_VCF.drop(columns=['QUAL', 'FILTER', 'INFO', 'FORMAT'], inplace=True)
subset_case_phased_VCF = phased_VCF[phased_VCF.ID.isin(list(CaseQV_pair.values()) + common_snps)]
# drop samples (columns) if there are only wild types in position where top significant locs
subset_phased_VCF = subset_case_phased_VCF[subset_case_phased_VCF.ID.isin(common_snps)]
columns2keep = subset_phased_VCF.loc[:, subset_phased_VCF.ne('0|0').any()].columns
subset_case_phased_VCF = subset_case_phased_VCF[columns2keep]
# drop snps (rows) if there are only wild types in position
subset_case_phased_VCF = subset_case_phased_VCF[subset_case_phased_VCF.iloc[:,5:].ne('0|0').any(axis=1)]

# drop samples have less than 3 snps
for sample_name in subset_case_phased_VCF.columns[~subset_case_phased_VCF.columns.isin(['CHROM', 'POS', 'ID', 'REF', 'ALT'])]:
    if subset_case_phased_VCF[sample_name].value_counts()['0|0'] > subset_case_phased_VCF.shape[0] - 3:
        subset_case_phased_VCF.drop(columns=sample_name, inplace=True)
    elif '1|0' not in subset_case_phased_VCF[sample_name].value_counts():
        subset_case_phased_VCF[sample_name].replace({'0|1': '1|0'}, inplace=True)
    elif '0|1' in subset_case_phased_VCF[sample_name].value_counts() and subset_case_phased_VCF[sample_name].value_counts()['0|1'] > subset_case_phased_VCF[sample_name].value_counts()['1|0']:
        subset_case_phased_VCF[sample_name].replace({'0|1': '1|0', '1|0': '0|1'}, inplace=True)


# internalSentinelSNPSampleName = list(subset_case_phased_VCF.columns[(subset_case_phased_VCF[subset_case_phased_VCF.ID == 'chr3-44804157-T-C'].isin(['1|1', '0|1', '1|0'])).iloc[0]])
# sentinelSNPSampleName = list(subset_case_phased_VCF.columns[(subset_case_phased_VCF[subset_case_phased_VCF.ID == 'chr3-44860894-T-A'].isin(['1|1', '0|1', '1|0'])).iloc[0]])
# with open('dbGAP/sampleData/internalSentinelSNPSampleName.txt', 'w') as filehandle:
#     for listitem in internalSentinelSNPSampleName:
#         filehandle.write('%s\n' % listitem)
# with open('dbGAP/sampleData/sentinelSNPSampleName.txt', 'w') as filehandle:
#     for listitem in sentinelSNPSampleName:
#         filehandle.write('%s\n' % listitem)

with open("dbGAP/sampleData/sentinelSNPSampleName.txt", "r") as file:
    sentinelSNPSampleName = file.readlines()
sentinelSNPSampleName = [x.strip() for x in sentinelSNPSampleName]

with open("dbGAP/sampleData/internalSentinelSNPSampleName.txt", "r") as file:
    internalSentinelSNPSampleName = file.readlines()
internalSentinelSNPSampleName = [x.strip() for x in internalSentinelSNPSampleName]

#QV samples with QV snps
QV_phased_vcf = subset_case_phased_VCF[subset_case_phased_VCF.ID.isin(CaseQV_pair.values())]
QV_phased_vcf = QV_phased_vcf[QV_phased_vcf.columns[~QV_phased_vcf.columns.isin(['CHROM', 'ID', 'REF', 'ALT'])]]
QV_phased_vcf = QV_phased_vcf.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="KIF15QVs")
QV_phased_vcf = QV_phased_vcf[QV_phased_vcf.KIF15QVs != "0|0"]

#drop QV samples
for sample in list(QV_phased_vcf['Sample Name']):
    if sample in sentinelSNPSampleName:
        sentinelSNPSampleName.remove(sample)
    if sample in internalSentinelSNPSampleName:
        internalSentinelSNPSampleName.remove(sample)

#QV samples with TopSNPs
topSNPs_phased_vcf = subset_case_phased_VCF[~subset_case_phased_VCF.ID.isin(CaseQV_pair.values())]
QV_samples = topSNPs_phased_vcf[['POS'] + list(QV_phased_vcf['Sample Name'])]
QV_samples = QV_samples.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
QV_samples = QV_samples[QV_samples.TopSNPs != "0|0"]

topSNPs_phased_vcf.drop(columns=list(QV_phased_vcf['Sample Name']) + ['CHROM', 'REF', 'ALT', 'ID'], inplace=True)
noneSentinelSNP = topSNPs_phased_vcf.drop(columns=list(np.union1d(internalSentinelSNPSampleName,sentinelSNPSampleName)))
noneSentinelSNP = noneSentinelSNP.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
noneSentinelSNP = noneSentinelSNP[noneSentinelSNP.TopSNPs != "0|0"]
# number of blue boxes
eitherSentinelSNP = topSNPs_phased_vcf[['POS'] + list(np.union1d(internalSentinelSNPSampleName,sentinelSNPSampleName))]
count_sum_df = pd.DataFrame()
count_sum_df['blue_count'] = (eitherSentinelSNP == '1|0').sum(axis=0)
# count_sum_df['orange_count'] = (eitherSentinelSNP == '0|1').sum(axis=0)
# count_sum_df['green_count'] = (eitherSentinelSNP == '1|1').sum(axis=0)
# count_sum_df['blue_orange_count'] = count_sum_df['blue_count'] + count_sum_df['orange_count']
ordered_columns = list(count_sum_df.sort_values(by = ['blue_count'], axis=0, ascending=True).index)
eitherSentinelSNP = eitherSentinelSNP[ordered_columns]
eitherSentinelSNP = eitherSentinelSNP.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
eitherSentinelSNP = eitherSentinelSNP[eitherSentinelSNP.TopSNPs != "0|0"]

# visualization
fig = plt.figure()
ax = fig.add_gridspec(33, 10)
ax1 = fig.add_subplot(ax[0, 0:10])
# pip install dna_features_viewer
from dna_features_viewer import GraphicFeature, GraphicRecord
features=[
    GraphicFeature(start=44804157, end=44804157, strand=+1, color="#ffd700",
                   label="rs74341405"),
    GraphicFeature(start=44860894, end=44860894, strand=+1, color="#ffd700",
                   label="rs78238620"),
    GraphicFeature(start=44761794, end=44853256, strand=+1, color="#ffcccc",
                   label="KIF15")]
record = GraphicRecord(sequence_length=160000, features=features, first_index=44710000)
record.plot(ax1)
ax1.margins(0)

ax1 = fig.add_subplot(ax[1:, 0:10])
ax1.set_xlim([44710000,44870000])
for dataset in noneSentinelSNP, eitherSentinelSNP, QV_samples:
    ax1.plot(dataset[dataset.TopSNPs == '1|1'].POS,
             dataset[dataset.TopSNPs == '1|1']['Sample Name'],
             mfc='lightgreen', marker='s', mec='forestgreen', ms=8, linestyle='none', label='KIF15 top SNPs - 1|1')
    ax1.plot(dataset[dataset.TopSNPs == '1|0'].POS,
            dataset[dataset.TopSNPs == '1|0']['Sample Name'],
            mfc='skyblue', marker='s', mec='blue', ms=8, linestyle='none', label='KIF15 top SNPs - 1|0')
    ax1.plot(dataset[dataset.TopSNPs == '0|1'].POS,
            dataset[dataset.TopSNPs == '0|1']['Sample Name'],
            mfc='orange', marker='s', mec='orangered', ms=8, linestyle='none', label='KIF15 top SNPs - 0|1')

ax1.plot(QV_phased_vcf[QV_phased_vcf.KIF15QVs == '1|0'].POS,
        QV_phased_vcf[QV_phased_vcf.KIF15QVs == '1|0']['Sample Name'],
        mfc='skyblue', marker='^', mec='blue', ms=8, linestyle='none', label='KIF15 QVs - 1|0')
ax1.plot(QV_phased_vcf[QV_phased_vcf.KIF15QVs == '0|1'].POS,
        QV_phased_vcf[QV_phased_vcf.KIF15QVs == '0|1']['Sample Name'],
        mfc='orange', marker='^', mec='orangered', ms=8, linestyle='none', label='KIF15 QVs - 0|1')

# customise
positions = np.arange(44720000, 44860001, 20000).tolist()
labels = ' '.join('{:,}'.format(x) for x in positions).split()
ax1.xaxis.set_major_locator(ticker.FixedLocator(positions))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.margins(0.01)
plt.grid(ls='--')
plt.subplots_adjust(left=0.125,
                    bottom=0.01,
                    right=0.9,
                    top=0.99,
                    wspace=0.2,
                    hspace=0.5)
fig.set_size_inches(28, 40)
fig.savefig('dbGAP_phasing_wSampleID.pdf', dpi=400)

# hide sampleID
ax1.set_yticklabels([])
fig.savefig('dbGAP_phasing_woSampleID.pdf', dpi=400)