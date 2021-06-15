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


with open('data/CaseSampleName_ID_dict_hg38.txt') as dict_file:
    CaseQV_pair = dict(x.rstrip().split(None, 1) for x in dict_file)

phased_VCF = read_vcf('dbGAP_IPF.KIF15.ID_nomono.vcf')
with open('data/KIF15_topsnps_hg38_longlist.txt', 'r') as f:
    common_snps = f.readlines()
# remove whitespace characters like `\n` at the end of each line
common_snps = [x.strip() for x in common_snps]

sentinel_snps = ['chr3-44804157-T-C', 'chr3-44860894-T-A']
phased_VCF[phased_VCF.ID.isin(list(CaseQV_pair.values()) + sentinel_snps)][['ID'] + list(CaseQV_pair.keys())].to_csv("KIF15_variants/KIF15_QV_sentinelSNPs_carrier.csv", index=False)

# [col for col in phased_VCF.columns if col in ['CHROM', 'POS', 'ID'] or "NWD" in col]
phased_VCF.drop(columns=['CHROM', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], inplace=True)
phased_VCF = phased_VCF[phased_VCF.ID.isin(list(CaseQV_pair.values()) + common_snps)]
# drop sample that are wild types
columns2keep = phased_VCF.loc[:, phased_VCF.ne('0|0').any()].columns
phased_VCF = phased_VCF[columns2keep]

# plot QV carriers only
KIF15QVCarrierVCF = phased_VCF[phased_VCF.ID.isin(list(CaseQV_pair.values()) + common_snps)][['POS', 'ID'] + list(CaseQV_pair.keys())]
for sample_name in KIF15QVCarrierVCF.columns[~KIF15QVCarrierVCF.columns.isin(['CHROM', 'POS', 'ID', 'REF', 'ALT'])]:
    if '1|0' not in KIF15QVCarrierVCF[sample_name].value_counts():
        KIF15QVCarrierVCF[sample_name].replace({'0|1': '1|0'}, inplace=True)
    elif '0|1' in KIF15QVCarrierVCF[sample_name].value_counts() and KIF15QVCarrierVCF[sample_name].value_counts()['0|1'] > KIF15QVCarrierVCF[sample_name].value_counts()['1|0']:
        KIF15QVCarrierVCF[sample_name].replace({'0|1': '1|0', '1|0': '0|1'}, inplace=True)

QVKIF15carrier = KIF15QVCarrierVCF[KIF15QVCarrierVCF.ID.isin(list(CaseQV_pair.values()))]
commonKIF15carrier = KIF15QVCarrierVCF[KIF15QVCarrierVCF.ID.isin(common_snps)]

def reformat2vis(input_vcf):
    input_vcf = input_vcf[[col for col in input_vcf.columns if col in ['POS'] or "NWD" in col]]
    input_vcf = input_vcf.melt(id_vars=["POS"],
                                 var_name="Sample Name",
                                 value_name="TopSNPs")
    input_vcf = input_vcf[input_vcf.TopSNPs != "0|0"]
    return(input_vcf)

QVKIF15carrier = reformat2vis(QVKIF15carrier)
commonKIF15carrier = reformat2vis(commonKIF15carrier)
# for sample_name in phased_VCF.columns[~phased_VCF.columns.isin(['CHROM', 'POS', 'ID', 'REF', 'ALT'])]:
#     if phased_VCF[sample_name].value_counts()['0|0'] > phased_VCF.shape[0] - 3:
#         phased_VCF.drop(columns=sample_name, inplace=True)
#     elif '1|0' not in phased_VCF[sample_name].value_counts():
#         phased_VCF[sample_name].replace({'0|1': '1|0'}, inplace=True)
#     elif '0|1' in phased_VCF[sample_name].value_counts() and phased_VCF[sample_name].value_counts()['0|1'] > phased_VCF[sample_name].value_counts()['1|0']:
#         phased_VCF[sample_name].replace({'0|1': '1|0', '1|0': '0|1'}, inplace=True)

fig = plt.figure()
ax = fig.add_gridspec(ncols=9, nrows=4)
ax1 = fig.add_subplot(ax[0, 0:9])
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
# ax1.xaxis.label.set_color('black')
# ax1.xaxis.tick_top()

ax1 = fig.add_subplot(ax[1:, 0:9])
ax1.set_xlim([44680000, 44880000])
ax1.plot(commonKIF15carrier[commonKIF15carrier.TopSNPs == '1|0'].POS,
        commonKIF15carrier[commonKIF15carrier.TopSNPs == '1|0']['Sample Name'],
        mfc='skyblue', marker='s', mec='blue', ms=8, linestyle='none', label='KIF15 top SNPs - 1|0')
ax1.plot(commonKIF15carrier[commonKIF15carrier.TopSNPs == '1|1'].POS,
        commonKIF15carrier[commonKIF15carrier.TopSNPs == '1|1']['Sample Name'],
        mfc='lightgreen', marker='s', mec='forestgreen', ms=8, linestyle='none', label='KIF15 top SNPs - 1|1')
ax1.plot(commonKIF15carrier[commonKIF15carrier.TopSNPs == '0|1'].POS,
        commonKIF15carrier[commonKIF15carrier.TopSNPs == '0|1']['Sample Name'],
        mfc='orange', marker='s', mec='orangered', ms=8, linestyle='none', label='KIF15 top SNPs - 0|1')

ax1.plot(QVKIF15carrier[QVKIF15carrier.TopSNPs == '1|0'].POS,
        QVKIF15carrier[QVKIF15carrier.TopSNPs == '1|0']['Sample Name'],
        mfc='skyblue', marker='^', mec='blue', ms=8, linestyle='none', label='KIF15 QVs - 1|0')
ax1.plot(QVKIF15carrier[QVKIF15carrier.TopSNPs == '0|1'].POS,
        QVKIF15carrier[QVKIF15carrier.TopSNPs == '0|1']['Sample Name'],
        mfc='orange', marker='^', mec='orangered', ms=8, linestyle='none', label='KIF15 QVs - 0|1')

positions = np.arange(44700000, 44880000, 25000).tolist()
labels = ' '.join('{:,}'.format(x) for x in positions).split()
ax1.xaxis.set_major_locator(ticker.FixedLocator(positions))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.margins(0.02)
plt.grid(ls='--')
ax1.xaxis.label.set_color('gray')

plt.subplots_adjust(left=0.125,
                    bottom=0.01,
                    right=0.9,
                    top=0.99,
                    wspace=0.2,
                    hspace=0.5)
fig.set_size_inches(18,8)
fig.savefig('KIF15QVCarrierVCF.pdf', dpi=400)

ax1.set_frame_on(False)
fig.set_size_inches(18,8)
fig.savefig('KIF15QVCarrierVCFwoFrame.pdf', dpi=400)