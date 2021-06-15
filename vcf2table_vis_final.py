import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import ticker

phased_vcf = pd.read_table("IPF_shortlist.KIF15.topsnps_longlist_AndQV.phased.txt")
with open('CaseSampleName_ID_dict.txt') as dict_file:
    CaseQV_pair = dict(x.rstrip().split(None, 1) for x in dict_file)

# only keep snps on the bigger blocks
phased_vcf = phased_vcf[phased_vcf.POS > 44756200]

for sample_name in phased_vcf.columns[~phased_vcf.columns.isin(['CHROM', 'POS', 'ID', 'REF', 'ALT'])]:
    if phased_vcf[sample_name].value_counts()['0|0'] > phased_vcf.shape[0] - 3:
        phased_vcf.drop(columns=sample_name, inplace=True)
    elif '1|0' not in phased_vcf[sample_name].value_counts():
        phased_vcf[sample_name].replace({'0|1': '1|0'}, inplace=True)
    elif '0|1' in phased_vcf[sample_name].value_counts() and phased_vcf[sample_name].value_counts()['0|1'] > phased_vcf[sample_name].value_counts()['1|0']:
        phased_vcf[sample_name].replace({'0|1': '1|0', '1|0': '0|1'}, inplace=True)


with open("sampleData/sentinelSNPSampleName.txt", "r") as file:
    sentinelSNPSampleName = file.readlines()
sentinelSNPSampleName = [x.strip() for x in sentinelSNPSampleName]

with open("sampleData/internalSentinelSNPSampleName.txt", "r") as file:
    internalSentinelSNPSampleName = file.readlines()
internalSentinelSNPSampleName = [x.strip() for x in internalSentinelSNPSampleName]

#QV samples with QV snps
QV_phased_vcf = phased_vcf[phased_vcf.ID.isin(CaseQV_pair.values())]
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
topSNPs_phased_vcf = phased_vcf[~phased_vcf.ID.isin(CaseQV_pair.values())]
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
count_sum_df['white_count'] = (eitherSentinelSNP == '0|0').sum(axis=0)
count_sum_df['orange_count'] = (eitherSentinelSNP == '0|1').sum(axis=0)
count_sum_df['green_count'] = (eitherSentinelSNP == '1|1').sum(axis=0)
# count_sum_df['blue_orange_count'] = count_sum_df['blue_count'] + count_sum_df['orange_count']
ordered_columns = list(count_sum_df.sort_values(by = ['green_count', 'orange_count', 'white_count'], axis=0, ascending=False).index)
eitherSentinelSNP = eitherSentinelSNP[ordered_columns]
eitherSentinelSNP = eitherSentinelSNP.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
eitherSentinelSNP = eitherSentinelSNP[eitherSentinelSNP.TopSNPs != "0|0"]

# visualization
fig = plt.figure()
ax = fig.add_gridspec(33, 10)
ax1 = fig.add_subplot(ax[0, 0:10])
from dna_features_viewer import GraphicFeature, GraphicRecord
features=[
    GraphicFeature(start=44845649, end=44845649, strand=+1, color="#ffd700",
                   label="rs74341405"),
    GraphicFeature(start=44902386, end=44902386, strand=+1, color="#ffd700",
                   label="rs78238620"),
    GraphicFeature(start=44803209, end=44894753, strand=+1, color="#ffcccc",
                   label="KIF15")]
record = GraphicRecord(sequence_length=160000, features=features, first_index=44750000)
record.plot(ax1)
ax1.margins(0)
# ax1.xaxis.label.set_color('black')
# ax1.xaxis.tick_top()

ax1 = fig.add_subplot(ax[1:, 0:10])
ax1.set_xlim([44750000,44910000])
for dataset in noneSentinelSNP, eitherSentinelSNP, QV_samples:
    ax1.plot(dataset[dataset.TopSNPs == '1|0'].POS,
            dataset[dataset.TopSNPs == '1|0']['Sample Name'],
            mfc='skyblue', marker='s', mec='blue', ms=8, linestyle='none', label='KIF15 top SNPs - 1|0')
    ax1.plot(dataset[dataset.TopSNPs == '1|1'].POS,
            dataset[dataset.TopSNPs == '1|1']['Sample Name'],
            mfc='lightgreen', marker='s', mec='forestgreen', ms=8, linestyle='none', label='KIF15 top SNPs - 1|1')
    ax1.plot(dataset[dataset.TopSNPs == '0|1'].POS,
            dataset[dataset.TopSNPs == '0|1']['Sample Name'],
            mfc='orange', marker='s', mec='orangered', ms=8, linestyle='none', label='KIF15 top SNPs - 0|1')

ax1.plot(QV_phased_vcf[QV_phased_vcf.KIF15QVs == '1|0'].POS,
        QV_phased_vcf[QV_phased_vcf.KIF15QVs == '1|0']['Sample Name'],
        mfc='skyblue', marker='^', mec='blue', ms=8, linestyle='none', label='KIF15 QVs - 1|0')
ax1.plot(QV_phased_vcf[QV_phased_vcf.KIF15QVs == '0|1'].POS,
        QV_phased_vcf[QV_phased_vcf.KIF15QVs == '0|1']['Sample Name'],
        mfc='orange', marker='^', mec='orangered', ms=8, linestyle='none', label='KIF15 QVs - 0|1')

positions = [44760000, 44780000, 44800000, 44820000, 44840000, 44860000, 44880000, 44900000]
labels = ['44,760,000', '44,780,000', '44,800,000', '44,820,000', '44,840,000', '44,860,000', '44,880,000', '44,900,000']
ax1.xaxis.set_major_locator(ticker.FixedLocator(positions))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.margins(0.01)
plt.grid(ls='--')
# ax1.xaxis.label.set_color('gray')

plt.subplots_adjust(left=0.125,
                    bottom=0.01,
                    right=0.9,
                    top=0.99,
                    wspace=0.2,
                    hspace=0.5)
# fig.legend(loc='center right')
fig.set_size_inches(28, 40)
fig.savefig('dbGAP_BigBlock_QV_final_dash.pdf', dpi=400)
