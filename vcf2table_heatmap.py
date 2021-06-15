import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

phased_vcf = pd.read_table("IPF_shortlist.KIF15.topsnps_longlist_AndQV.phased.txt")

vcf_df = pd.Series()
for sample_name in phased_vcf.columns[5:]:
    strandA = phased_vcf[sample_name].str.split("|", expand=True).astype(int)[0]
    strandB = phased_vcf[sample_name].str.split("|", expand=True).astype(int)[1]
    vcf_df = vcf_df.append(strandA[strandA == 1])
    vcf_df[vcf_df == 1] = sample_name + "_A"
    vcf_df = vcf_df.append(strandB[strandB == 1])
    vcf_df[vcf_df == 1] = sample_name + "_B"

import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, plt = plt.subplots()
ax1.set_xlim([44700000,44910000])
ax1.get_xax1is().tick_bottom()
# ax1.ax1es.get_yax1is().set_visible(False)
ax1.set_frame_on(False)
ax1.grid(color='gray', linestyle='-', linewidth=1, ax1is='y')
ax1.plot(vcf_df, mfc='navy', marker='o', mec='blue', ms=8, linestyle='none', label='KIF15 top SNPs')
fig.set_size_inches(20, 90)
fig.savefig('test_longIPFlist.pdf', dpi=100)

## second way
phased_vcf = pd.read_table("IPF_shortlist.KIF15.topsnps_longlist_AndQV.phased.txt")
topSNPs_phased_vcf = phased_vcf[~phased_vcf.ID.isin(CaseQV_pair.values())]
QV_phased_vcf = phased_vcf[phased_vcf.ID.isin(CaseQV_pair.values())]
topSNPs_phased_vcf = topSNPs_phased_vcf.iloc[:, [1] + list(range(5,441))]
topSNPs_phased_vcf = topSNPs_phased_vcf.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="TopSNPs")
topSNPs_phased_vcf = topSNPs_phased_vcf[topSNPs_phased_vcf.TopSNPs != "0|0"]
QV_phased_vcf = QV_phased_vcf.iloc[:, [1] + list(range(5,441))]
QV_phased_vcf = QV_phased_vcf.melt(id_vars=["POS"],
        var_name="Sample Name",
        value_name="KIF15QVs")
QV_phased_vcf = QV_phased_vcf[QV_phased_vcf.KIF15QVs != "0|0"]

fig, (ax1, ax2) = plt.subplots(2)
ax1.set_xlim([44720000,44910000])
ax1.figure.set_size_inches(25, 80)
# ax1.get_xax1is().tick_bottom()
# ax1.ax1es.get_yax1is().set_visible(False)
# ax1.set_frame_on(False)
# sns.scatterplot(x='POS', y='Sample Name', data=topSNPs_phased_vcf, hue='TopSNPs', ec=None, ax1=ax1)
# sns.scatterplot(x='POS', y='Sample Name', data=QV_phased_vcf, hue='KIF15QVs', ec=None, ax1=ax1)
ax1.plot(topSNPs_phased_vcf[topSNPs_phased_vcf.TopSNPs == '1|0'].POS,
        topSNPs_phased_vcf[topSNPs_phased_vcf.TopSNPs == '1|0']['Sample Name'],
        mfc='skyblue', marker='s', mec='blue', ms=8, linestyle='none', label='KIF15 top SNPs - 1|0')
ax1.plot(topSNPs_phased_vcf[topSNPs_phased_vcf.TopSNPs == '0|1'].POS,
        topSNPs_phased_vcf[topSNPs_phased_vcf.TopSNPs == '0|1']['Sample Name'],
        mfc='orange', marker='s', mec='orangered', ms=8, linestyle='none', label='KIF15 top SNPs - 0|1')
ax1.plot(topSNPs_phased_vcf[topSNPs_phased_vcf.TopSNPs == '1|1'].POS,
        topSNPs_phased_vcf[topSNPs_phased_vcf.TopSNPs == '1|1']['Sample Name'],
        mfc='lightgreen', marker='s', mec='forestgreen', ms=8, linestyle='none', label='KIF15 top SNPs - 1|1')

ax1.plot(QV_phased_vcf[QV_phased_vcf.KIF15QVs == '1|0'].POS,
        QV_phased_vcf[QV_phased_vcf.KIF15QVs == '1|0']['Sample Name'],
        mfc='skyblue', marker='^', mec='blue', ms=8, linestyle='none', label='KIF15 QVs - 1|0')
ax1.plot(QV_phased_vcf[QV_phased_vcf.KIF15QVs == '0|1'].POS,
        QV_phased_vcf[QV_phased_vcf.KIF15QVs == '0|1']['Sample Name'],
        mfc='orange', marker='^', mec='orangered', ms=8, linestyle='none', label='KIF15 QVs - 0|1')


from dna_features_viewer import GraphicFeature, GraphicRecord
features=[
    GraphicFeature(start=44902386, end=44902386, strand=+1, color="#ffd700",
                   label="rs78238620"),
    GraphicFeature(start=44803209, end=44894753, strand=+1, color="#ffcccc",
                   label="KIF15")]
record = GraphicRecord(sequence_length=190000, features=features, first_index=44720000)
record.plot(ax2)

plt.subplots_adjust(left=0.125,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.2,
                    hspace=0)
fig.legend(loc='center right')
fig.savefig('test_QV_longIPFlist.pdf', dpi=100)