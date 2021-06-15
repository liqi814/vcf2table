import pandas as pd
import io
import numpy as np

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


with open('CaseSampleName_ID_dict.txt') as dict_file:
    CaseQV_pair = dict(x.rstrip().split(None, 1) for x in dict_file)

phased_VCF = read_vcf('IPF.KIF15.500kb.3410s.ID.bi-allelic_only_nomiss.recode.phased.vcf')
with open('KIF15_topsnps_hg19_longlist.txt', 'r') as f:
    common_snps = f.readlines()
# remove whitespace characters like `\n` at the end of each line
common_snps = [x.strip() for x in common_snps]
with open('Gundula_after_pruning_All_case_1c.txt', 'r') as f:
    Case_name = f.readlines()
# remove whitespace characters like `\n` at the end of each line
# palm samples are WES data
Case_name = [x.strip() for x in Case_name if 'palm' not in x]

subset_case_phased_VCF = phased_VCF[['CHROM', 'POS', 'ID', 'REF', 'ALT'] + Case_name]
subset_case_phased_VCF = subset_case_phased_VCF[subset_case_phased_VCF.ID.isin(list(CaseQV_pair.values()) + common_snps)]
subset_case_phased_VCF.to_csv("IPF.KIF15.topsnps_longlist_AndQV.phased.txt", index=False, sep='\t')

#drop samples (columns) if there are only wild types in position where top significant locs
subset_phased_VCF = subset_case_phased_VCF[subset_case_phased_VCF.ID.isin(common_snps)]
columns2keep = subset_phased_VCF.loc[:, subset_phased_VCF.ne('0|0').any()].columns
subset_case_phased_VCF = subset_case_phased_VCF[columns2keep]
#drop snps (rows) if there are only wild types in position
subset_case_phased_VCF = subset_case_phased_VCF[subset_case_phased_VCF.iloc[:,5:].ne('0|0').any(ax1is=1)]
subset_case_phased_VCF.to_csv("IPF_shortlist.KIF15.topsnps_longlist_AndQV.phased.txt", index=False, sep='\t')

subset_case_phased_VCF.reset_index(drop=True, inplace=True)

column_names = ["POS", "ID", "StrandA", "StrandB", "StrandApos", "StrandBpos"]
sample_df = pd.DataFrame(columns = column_names)

for index, sample_name in enumerate(subset_case_phased_VCF.columns[5:]):
    sample_df2 = pd.DataFrame(
        {
            'POS': subset_case_phased_VCF['POS'].astype(int),
            'ID': subset_case_phased_VCF['ID'],
            'StrandA': subset_case_phased_VCF[sample_name].str.split("|", expand=True).astype(int)[0],
            'StrandB': subset_case_phased_VCF[sample_name].str.split("|", expand=True).astype(int)[1],
            'StrandApos': 0,
            'StrandBpos': 0
        })
    if sample_df2.StrandA.sum() < sample_df2.StrandB.sum():
        sample_df2.rename(columns={'StrandA': 'StrandB', 'StrandB': 'StrandA'}, inplace=True)
        sample_df2 = sample_df2[["POS", "ID", "StrandA", "StrandB", "StrandApos", "StrandBpos"]]
    sample_df = sample_df.append(sample_df2, ignore_index=True)

for variants_ID in sample_df['ID'].unique():
    StrandAposidx = 1.0
    StrandBposidx = -1.0
    for idx, rows in sample_df[sample_df.ID == variants_ID].iterrows():
        if rows['StrandA'] == 1:
            StrandAposidx += 0.5
            sample_df.loc[idx, 'StrandApos'] = StrandAposidx
        if rows['StrandB'] == 1:
            StrandBposidx -= 0.5
            sample_df.loc[idx, 'StrandBpos'] = StrandBposidx


# visualization
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import mpl_toolkits.ax1isartist as ax1isartist


x_limit = list(subset_case_phased_VCF.iloc[[0,-1]]['POS'])

KIF15_coordinates = [44803209,44894753]

# for index, row in subset_phased_VCF[['POS','IPF0090']].iterrows():
#     print(row['POS'])
# fig, ax1 = plt.subplots()
# ax1.set_frame_on(False)

fig = plt.figure()
ax1 = ax1isartist.Subplot(fig, 111)
fig.add_ax1es(ax1)
ax1.set_xlim([44700000,44910000])
# ax1.set_frame_on(False)
# remove y ticks
ax1.ax1is["bottom"].set_ax1isline_style("-|>", size = 1.5)
ax1.ax1is["left"].set_ax1isline_style("->", size = 1.5)
plt.yticks([])
# plt.xticks([])
# #通过set_visible方法设置绘图区的顶部及右侧坐标轴隐藏
ax1.ax1is["top"].set_visible(False)
ax1.ax1is["right"].set_visible(False)

ax1.hlines(y=1, xmin=44700000, xmax1=44910000, linewidth=2, color='gray', linestyles='dashed')
ax1.hlines(y=-1, xmin=44700000, xmax1=44910000, linewidth=2, color='gray', linestyles='dashed')

ax1.plot(sample_df[sample_df.StrandA == 1]['POS'], sample_df[sample_df.StrandA == 1]['StrandApos'],
        mfc='lightblue', marker='o', mec='lightblue', ms=2, linestyle='none', label='KIF15 top SNPs')
ax1.plot(sample_df[sample_df.StrandB == 1]['POS'], sample_df[sample_df.StrandB == 1]['StrandBpos'],
        mfc='lightblue', marker='o', mec='lightblue', ms=2, linestyle='none', label='KIF15 top SNPs')

plt.show()
ax1.hlines(y=-1, xmin=44800000, xmax1=44905000, linewidth=2, color='black')
# ax1.add_patch(patches.Rectangle((44803209, 0), 91544, 0.3, color='none', ec='black', label='KIF15'))
ax1.add_patch(patches.Rectangle((44803209, -1), 91544, 1, color='none', ec='black'))
ax1.text(np.mean(KIF15_coordinates), -0.5, 'KIF15', style='italic', ha='center', va='center')
ax1.plot(44902386, -1, mfc='black', marker='o', mec='black', ms=5)
ax1.text(44902386, -0.5, 'rs78238620', va='center')

column_names = ["POS", "ID", "LEFT", "RIGHT", "LEFTidx", "RIGHTidx"]
sample_df = pd.DataFrame(columns = column_names)
QV_Carrier_df = pd.DataFrame(columns = column_names)

for index, sample_name in enumerate(subset_phased_VCF.columns[5:]):
    ax1.hlines(y=3*index+1, xmin=44800000, xmax1=44905000, linewidth=2, color='black')
    ax1.hlines(y=3*index+2, xmin=44800000, xmax1=44905000, linewidth=2, color='black')
    ax1.text(44800000, 3 * index + 1.5, sample_name, ha='right', va='center')
    sample_df2 = pd.DataFrame(
        {
            'POS': subset_phased_VCF['POS'].astype(int),
            'ID': subset_phased_VCF['ID'],
            'LEFT': subset_phased_VCF[sample_name].str.split("|", expand=True).astype(int)[0],
            'RIGHT': subset_phased_VCF[sample_name].str.split("|", expand=True).astype(int)[1],
            'LEFTidx': [3 * index + 1] * subset_phased_VCF.shape[0],
            'RIGHTidx': [3 * index + 2] * subset_phased_VCF.shape[0]
        }
    )
    sample_df = sample_df.append(sample_df2, ignore_index=True)

# find KIF15 QV snps for each QV carrier
QV_Carrier = sample_df[sample_df['ID'].isin(list(CaseQV_pair.values()))].loc[(sample_df[['LEFT', 'RIGHT']]==1).any(ax1is=1), :]
# keep all KIF15 topsnps in a new variable
KIF15_topsnps = sample_df[~sample_df['ID'].isin(list(CaseQV_pair.values()))]
# reshape data frame
QV_Carrier = QV_Carrier[['POS', 'ID', 'LEFT', 'LEFTidx']].rename(columns = {'LEFT':'Alleles', 'LEFTidx':'POSidx'}).append(
    QV_Carrier[['POS', 'ID', 'RIGHT', 'RIGHTidx']].rename(columns = {'RIGHT':'Alleles', 'RIGHTidx':'POSidx'}), ignore_index=True)
QV_Carrier = QV_Carrier.sort_values('POS')
KIF15_topsnps = KIF15_topsnps[['POS', 'ID', 'LEFT', 'LEFTidx']].rename(columns = {'LEFT':'Alleles', 'LEFTidx':'POSidx'}).append(
    KIF15_topsnps[['POS', 'ID', 'RIGHT', 'RIGHTidx']].rename(columns = {'RIGHT':'Alleles', 'RIGHTidx':'POSidx'}), ignore_index=True)
KIF15_topsnps = KIF15_topsnps.sort_values('POS')

# ax1.plot(KIF15_topsnps[KIF15_topsnps.Alleles == 0]['POS'], KIF15_topsnps[KIF15_topsnps.Alleles == 0]['POSidx'],
#         mfc='white', marker='o', mec='blue', ms=8, linestyle='none', label='Wild type')
# ax1.plot(QV_Carrier[QV_Carrier.Alleles == 0]['POS'], QV_Carrier[QV_Carrier.Alleles == 0]['POSidx'],
#         mfc='white', marker='o', mec='orangered', ms=8, linestyle='none', label='Wild type')
ax1.plot(KIF15_topsnps[KIF15_topsnps.Alleles == 1]['POS'], KIF15_topsnps[KIF15_topsnps.Alleles == 1]['POSidx'],
        mfc='navy', marker='o', mec='blue', ms=8, linestyle='none', label='KIF15 top SNPs')
ax1.plot(QV_Carrier[QV_Carrier.Alleles == 1]['POS'], QV_Carrier[QV_Carrier.Alleles == 1]['POSidx'],
        mfc='orange', marker='o', mec='orangered', ms=8, linestyle='none', label='KIF15 QVs')
fig.legend(loc='center right')
fig.set_size_inches(20, 16)
fig.savefig('test2pngNoWildTypelowdpi.pdf', dpi=100)