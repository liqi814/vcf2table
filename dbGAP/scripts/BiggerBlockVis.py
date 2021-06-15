from scripts.vcf2table import CaseQV_pair, common_snps, phased_VCF, reformat2vis
import numpy as np

bigger_phased_vcf = phased_VCF[phased_VCF.POS > 44710000]

for sample_name in bigger_phased_vcf.columns[~bigger_phased_vcf.columns.isin(['CHROM', 'POS', 'ID', 'REF', 'ALT'])]:
    if bigger_phased_vcf[sample_name].value_counts()['0|0'] > bigger_phased_vcf.shape[0] - 3:
        bigger_phased_vcf.drop(columns=sample_name, inplace=True)
    elif '1|0' not in bigger_phased_vcf[sample_name].value_counts():
        bigger_phased_vcf[sample_name].replace({'0|1': '1|0'}, inplace=True)
    elif '0|1' in bigger_phased_vcf[sample_name].value_counts() and bigger_phased_vcf[sample_name].value_counts()['0|1'] > bigger_phased_vcf[sample_name].value_counts()['1|0']:
        bigger_phased_vcf[sample_name].replace({'0|1': '1|0', '1|0': '0|1'}, inplace=True)


sentinelSNPSampleName = bigger_phased_vcf.columns[(bigger_phased_vcf[bigger_phased_vcf.ID == 'chr3-44860894-T-A'] != '0|0').iloc[0]]
sentinelSNPSampleName = sentinelSNPSampleName.to_series()[~sentinelSNPSampleName.to_series().isin(['CHROM', 'POS', 'ID', 'REF', 'ALT'])]
sentinelSNPSampleName.to_csv("sampleData/sentinelSNPSampleName.txt", index=False, header=False)

internalSentinelSNPSampleName = bigger_phased_vcf.columns[(bigger_phased_vcf[bigger_phased_vcf.ID == 'chr3-44804157-T-C'] != '0|0').iloc[0]]
internalSentinelSNPSampleName = internalSentinelSNPSampleName.to_series()[~internalSentinelSNPSampleName.to_series().isin(['CHROM', 'POS', 'ID', 'REF', 'ALT'])]
internalSentinelSNPSampleName.to_csv("sampleData/internalSentinelSNPSampleName.txt", index=False, header=False)

morethan3SNPsSampleName = bigger_phased_vcf.columns[~bigger_phased_vcf.columns.isin(['CHROM', 'POS', 'ID', 'REF', 'ALT'])]
morethan3SNPsSampleName.to_series().to_csv("sampleData/morethan3SNPsSampleName.txt", index=False, header=False)

QV_phased_vcf = bigger_phased_vcf[bigger_phased_vcf.ID.isin(CaseQV_pair.values())]
QV_phased_vcf = QV_phased_vcf[QV_phased_vcf.columns[~QV_phased_vcf.columns.isin(['CHROM', 'ID', 'REF', 'ALT'])]]
QV_phased_vcf = QV_phased_vcf.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="KIF15QVs")
QV_phased_vcf = QV_phased_vcf[QV_phased_vcf.KIF15QVs != "0|0"]

internalSentinelSNPSampleName = internalSentinelSNPSampleName.to_numpy()
sentinelSNPSampleName = sentinelSNPSampleName.to_numpy()
morethan3SNPsSampleName = morethan3SNPsSampleName.to_numpy()
#drop QV samples
for sample in list(QV_phased_vcf['Sample Name']):
    if sample in sentinelSNPSampleName:
        sentinelSNPSampleName.remove(sample)
    if sample in internalSentinelSNPSampleName:
        internalSentinelSNPSampleName.remove(sample)

#QV samples
topSNPs_phased_vcf = bigger_phased_vcf[~bigger_phased_vcf.ID.isin(CaseQV_pair.values())]
QV_samples = topSNPs_phased_vcf[['POS'] + list(QV_phased_vcf['Sample Name'])]
QV_samples = QV_samples.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
QV_samples = QV_samples[QV_samples.TopSNPs != "0|0"]
topSNPs_phased_vcf.drop(columns=list(QV_phased_vcf['Sample Name']), inplace=True)
# fist plot the sample have both internal and external sentinel snps and reordered by internal snp and ['3-44902386-T-A', '3-44845649-T-C', '3-44786946-T-A']
# have both two snps
bothSentinelSNP = topSNPs_phased_vcf[['POS', 'ID'] + list(np.intersect1d(internalSentinelSNPSampleName,sentinelSNPSampleName))]
# reorder sample names for better visualization
bothSentinelSNP.set_index('ID', inplace=True)
tmp_phased_vcf = bothSentinelSNP[bothSentinelSNP.columns[~bothSentinelSNP.columns.isin(['CHROM', 'POS', 'REF', 'ALT'])]]
tmp_phased_vcf.sort_values(by = ['chr3-44804157-T-C', 'chr3-44745454-T-A'], axis = 1, ascending = False, inplace=True)
bothSentinelSNP.reset_index(inplace=True)
bothSentinelSNP = bothSentinelSNP[['POS'] + list(tmp_phased_vcf.columns)]
bothSentinelSNP = bothSentinelSNP.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
bothSentinelSNP = bothSentinelSNP[bothSentinelSNP.TopSNPs != "0|0"]
# doesn't have either internal and external snps
noneSentinelSNP = topSNPs_phased_vcf.drop(columns=['ID'] + list(np.union1d(internalSentinelSNPSampleName, sentinelSNPSampleName)))
noneSentinelSNP = noneSentinelSNP.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
noneSentinelSNP = noneSentinelSNP[noneSentinelSNP.TopSNPs != "0|0"]
# # has internal but not external
# internalSentinelSNP = topSNPs_phased_vcf[["POS"] + list(np.setdiff1d(internalSentinelSNPSampleName,sentinelSNPSampleName))]
# internalSentinelSNP = internalSentinelSNP.melt(id_vars=["POS"],
#         var_name="Sample Name",
#         value_name="TopSNPs")
# internalSentinelSNP = internalSentinelSNP[internalSentinelSNP.TopSNPs != "0|0"]
# has external but not internal
SentinelSNP = topSNPs_phased_vcf[["POS"] + list(np.setdiff1d(sentinelSNPSampleName, internalSentinelSNPSampleName))]
SentinelSNP = SentinelSNP.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
SentinelSNP = SentinelSNP[SentinelSNP.TopSNPs != "0|0"]

# visualization
fig = plt.figure()
ax = fig.add_gridspec(33, 10)
ax1 = fig.add_subplot(ax[0, 0:10])
from dna_features_viewer import GraphicFeature, GraphicRecord
features=[
    GraphicFeature(start=44804157, end=44804157, strand=+1, color="#ffd700",
                   label="rs74341405"),
    GraphicFeature(start=44860894, end=44860894, strand=+1, color="#ffd700",
                   label="rs78238620"),
    GraphicFeature(start=44761794, end=44853256, strand=+1, color="#ffcccc",
                   label="KIF15")]
record = GraphicRecord(sequence_length=165000, features=features, first_index=44710000)
record.plot(ax1)
ax1.margins(0)
# ax1.xaxis.label.set_color('black')
# ax1.xaxis.tick_top()

ax1 = fig.add_subplot(ax[1:, 0:10])
ax1.set_xlim([44710000,44875000])
for dataset in bothSentinelSNP, SentinelSNP, noneSentinelSNP, QV_samples:
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

positions = np.arange(44720000, 44860001, 20000).tolist()
labels = ' '.join('{:,}'.format(x) for x in positions).split()
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
fig.set_size_inches(30, 40)
fig.savefig('internalsentinelSNP_BigBlock_QV_longIPFlist_dash.pdf', dpi=400)