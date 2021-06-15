import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FuncFormatter
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

# # reorder sample names for better visualization
# phased_vcf.set_index('ID', inplace=True)
# tmp_phased_vcf = phased_vcf[phased_vcf.columns[~phased_vcf.columns.isin(['CHROM', 'POS', 'REF', 'ALT'])]]
# # tmp_phased_vcf.sort_values(by = ['3-44902386-T-A', '3-44845649-T-C'], axis = 1, ascending = False, inplace=True)
# tmp_phased_vcf.sort_values(by = ['3-44902386-T-A', '3-44845649-T-C', '3-44786946-T-A'], axis = 1, ascending = False, inplace=True)
# # tmp_phased_vcf.sort_values(by = '3-44902386-T-A', axis = 1, ascending = True, inplace=True)
# phased_vcf.reset_index(inplace=True)
# phased_vcf = phased_vcf[['CHROM', 'POS', 'REF','ID'] + list(tmp_phased_vcf.columns)]

with open("sampleData/sentinelSNPSampleName.txt", "r") as file:
    sentinelSNPSampleName = file.readlines()
sentinelSNPSampleName = [x.strip() for x in sentinelSNPSampleName]

with open("sampleData/internalSentinelSNPSampleName.txt", "r") as file:
    internalSentinelSNPSampleName = file.readlines()
internalSentinelSNPSampleName = [x.strip() for x in internalSentinelSNPSampleName]

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

#QV samples
topSNPs_phased_vcf = phased_vcf[~phased_vcf.ID.isin(CaseQV_pair.values())]
QV_samples = topSNPs_phased_vcf[['POS'] + list(QV_phased_vcf['Sample Name'])]
QV_samples = QV_samples.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
QV_samples = QV_samples[QV_samples.TopSNPs != "0|0"]
topSNPs_phased_vcf.drop(columns=list(QV_phased_vcf['Sample Name']), inplace=True)
# fist plot the sample have both internal and external sentinel snps and reordered by internal snp and ['3-44902386-T-A', '3-44845649-T-C', '3-44786946-T-A']
# have both two snps
bothSentinelSNP = topSNPs_phased_vcf[['CHROM', 'POS', 'ID', 'REF', 'ALT'] + list(np.intersect1d(internalSentinelSNPSampleName,sentinelSNPSampleName))]
# reorder sample names for better visualization
bothSentinelSNP.set_index('ID', inplace=True)
tmp_phased_vcf = bothSentinelSNP[bothSentinelSNP.columns[~bothSentinelSNP.columns.isin(['CHROM', 'POS', 'REF', 'ALT'])]]
tmp_phased_vcf.sort_values(by = ['3-44845649-T-C', '3-44786946-T-A'], axis = 1, ascending = False, inplace=True)
bothSentinelSNP.reset_index(inplace=True)
bothSentinelSNP = bothSentinelSNP[['POS'] + list(tmp_phased_vcf.columns)]
bothSentinelSNP = bothSentinelSNP.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
bothSentinelSNP = bothSentinelSNP[bothSentinelSNP.TopSNPs != "0|0"]
# doesn't have either internal and external snps
noneSentinelSNP = topSNPs_phased_vcf.drop(columns=['CHROM', 'ID', 'REF', 'ALT'] + list(np.union1d(internalSentinelSNPSampleName,sentinelSNPSampleName)))
noneSentinelSNP = noneSentinelSNP.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
noneSentinelSNP = noneSentinelSNP[noneSentinelSNP.TopSNPs != "0|0"]
# has internal but not external
internalSentinelSNP = topSNPs_phased_vcf[["POS"] + list(np.setdiff1d(internalSentinelSNPSampleName,sentinelSNPSampleName))]
internalSentinelSNP = internalSentinelSNP.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
internalSentinelSNP = internalSentinelSNP[internalSentinelSNP.TopSNPs != "0|0"]
# has external but not internal
SentinelSNP = topSNPs_phased_vcf[["POS"] + list(np.setdiff1d(sentinelSNPSampleName, internalSentinelSNPSampleName))]
SentinelSNP = SentinelSNP.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
SentinelSNP = SentinelSNP[SentinelSNP.TopSNPs != "0|0"]

# # put QV samples in the end and show on the top of figure
# topSNPs_phased_vcf = phased_vcf[~phased_vcf.ID.isin(CaseQV_pair.values())]
# for top_sample in list(QV_phased_vcf["Sample Name"]):
#     topSNPs_phased_vcf.insert(len(topSNPs_phased_vcf.columns)-1, top_sample, topSNPs_phased_vcf.pop(top_sample))
#
# topSNPs_phased_vcf = topSNPs_phased_vcf[topSNPs_phased_vcf.columns[~topSNPs_phased_vcf.columns.isin(['CHROM', 'ID', 'REF', 'ALT'])]]
# topSNPs_phased_vcf = topSNPs_phased_vcf.melt(id_vars=["POS"],
#         var_name="Sample Name",
#         value_name="TopSNPs")
# topSNPs_phased_vcf = topSNPs_phased_vcf[topSNPs_phased_vcf.TopSNPs != "0|0"]

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
for dataset in bothSentinelSNP, SentinelSNP, internalSentinelSNP, noneSentinelSNP, QV_samples:
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
fig.savefig('internalsentinelSNP_BigBlock_QV_longIPFlist_dash.pdf', dpi=400)
