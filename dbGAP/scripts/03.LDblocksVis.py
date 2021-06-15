# both LD blcoks
from scripts.vcf2table import CaseQV_pair, common_snps, phased_VCF, reformat2vis

for sample_name in phased_VCF.columns[~phased_VCF.columns.isin(['CHROM', 'POS', 'ID', 'REF', 'ALT'])]:
    if phased_VCF[sample_name].value_counts()['0|0'] > phased_VCF.shape[0] - 3:
        phased_VCF.drop(columns=sample_name, inplace=True)
    elif '1|0' not in phased_VCF[sample_name].value_counts():
        phased_VCF[sample_name].replace({'0|1': '1|0'}, inplace=True)
    elif '0|1' in phased_VCF[sample_name].value_counts() and phased_VCF[sample_name].value_counts()['0|1'] > phased_VCF[sample_name].value_counts()['1|0']:
        phased_VCF[sample_name].replace({'0|1': '1|0', '1|0': '0|1'}, inplace=True)

QVcarrier = phased_VCF[phased_VCF.ID.isin(list(CaseQV_pair.values()))]
commoncarrier = phased_VCF[phased_VCF.ID.isin(common_snps)]
commoncarrier.set_index('ID', inplace=True)
tmp_commoncarrier = commoncarrier[commoncarrier.columns[~commoncarrier.columns.isin(['CHROM', 'POS', 'REF', 'ALT'])]]
tmp_commoncarrier.sort_values(by = ['chr3-44860894-T-A', 'chr3-44804157-T-C', 'chr3-44685971-A-G'], axis = 1, ascending = False, inplace=True)
# tmp_commoncarrier.sort_values(by = '3-44902386-T-A', axis = 1, ascending = True, inplace=True)
commoncarrier.reset_index(inplace=True)
commoncarrier = commoncarrier[['POS','ID'] + list(tmp_commoncarrier.columns)]


QVcarrier = reformat2vis(QVcarrier)
commoncarrier = reformat2vis(commoncarrier)

# visualization
fig = plt.figure()
ax = fig.add_gridspec(33, 8)
ax1 = fig.add_subplot(ax[0, 0:8])
from dna_features_viewer import GraphicFeature, GraphicRecord
features=[
    GraphicFeature(start=44804157, end=44804157, strand=+1, color="#ffd700",
                   label="rs74341405"),
    GraphicFeature(start=44860894, end=44860894, strand=+1, color="#ffd700",
                   label="rs78238620"),
    GraphicFeature(start=44761794, end=44853256, strand=+1, color="#ffcccc",
                   label="KIF15")]
record = GraphicRecord(sequence_length=200000, features=features, first_index=44680000)
record.plot(ax1)
ax1.margins(0)

ax1 = fig.add_subplot(ax[1:, 0:8])
ax1.set_xlim([44680000, 44880000])
ax1.plot(commoncarrier[commoncarrier.TopSNPs == '1|0'].POS,
        commoncarrier[commoncarrier.TopSNPs == '1|0']['Sample Name'],
        mfc='skyblue', marker='s', mec='blue', ms=8, linestyle='none', label='KIF15 top SNPs - 1|0')
ax1.plot(commoncarrier[commoncarrier.TopSNPs == '0|1'].POS,
        commoncarrier[commoncarrier.TopSNPs == '0|1']['Sample Name'],
        mfc='orange', marker='s', mec='orangered', ms=8, linestyle='none', label='KIF15 top SNPs - 0|1')
ax1.plot(commoncarrier[commoncarrier.TopSNPs == '1|1'].POS,
        commoncarrier[commoncarrier.TopSNPs == '1|1']['Sample Name'],
        mfc='lightgreen', marker='s', mec='forestgreen', ms=8, linestyle='none', label='KIF15 top SNPs - 1|1')

ax1.plot(QVcarrier[QVcarrier.TopSNPs == '1|0'].POS,
        QVcarrier[QVcarrier.TopSNPs == '1|0']['Sample Name'],
        mfc='skyblue', marker='^', mec='blue', ms=8, linestyle='none', label='KIF15 QVs - 1|0')
ax1.plot(QVcarrier[QVcarrier.TopSNPs == '0|1'].POS,
        QVcarrier[QVcarrier.TopSNPs == '0|1']['Sample Name'],
        mfc='orange', marker='^', mec='orangered', ms=8, linestyle='none', label='KIF15 QVs - 0|1')

positions = np.arange(44700000, 44880000, 25000).tolist()
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
fig.set_size_inches(26, 40)
fig.savefig('BothBlocksLongIPFlist_dbGAP.pdf', dpi=400)