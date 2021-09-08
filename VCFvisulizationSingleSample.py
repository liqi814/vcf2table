import os
os.chdir("D:/All_IPF_CNV_analysis/candiate_CNV")
import io
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
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


vcf_file = read_vcf('IPF3170.197154_PARN_output_biallele.phased.vcf')
vcf_file = vcf_file[(vcf_file['IPF3170.197154'] != '0|0') & (vcf_file['POS'] >= 14440000) & (vcf_file['POS'] <= 14740000)].reset_index(drop=True)
vcf_file['Sample Name'] = pd.Series(["IPF3170.197154" for x in range(len(vcf_file.index))])

missenseQV = vcf_file[vcf_file.POS == 14647972]
otherSNPs = vcf_file[vcf_file.POS != 14647972]

# visualization
fig = plt.figure()
ax = fig.add_gridspec(3, 10)
ax1 = fig.add_subplot(ax[0, 0:10])
# pip install dna_features_viewer
from dna_features_viewer import GraphicFeature, GraphicRecord
features=[
    GraphicFeature(start=14647972, end=14647972, strand=+1, color="#ffd700",
                   label="c.1172T>C (p.Phe391Ser)"),
    GraphicFeature(start=14672001, end=14737000, strand=+1, color="#ffd700",
                   label="Candidate CNV"),
    GraphicFeature(start=14500000, end=14724128, strand=+1, color="#ffcccc",
                   label="PARN")]
record = GraphicRecord(sequence_length=240000, features=features, first_index=14500000)
record.plot(ax1)
ax1.margins(0)

ax1 = fig.add_subplot(ax[1:, 0:10])
ax1.set_xlim([14500000,14740000])
ax1.plot(otherSNPs[otherSNPs['IPF3170.197154'] == '1|1'].POS, otherSNPs[otherSNPs['IPF3170.197154'] == '1|1']['Sample Name'],
        mfc='skyblue', marker='s', mec='blue', ms=8, linestyle='none', label='1|0')
ax1.plot(otherSNPs[otherSNPs['IPF3170.197154'] == '1|0'].POS, otherSNPs[otherSNPs['IPF3170.197154'] == '1|0']['Sample Name'],
        mfc='lightgreen', marker='s', mec='forestgreen', ms=8, linestyle='none', label='1|1')
ax1.plot(otherSNPs[otherSNPs['IPF3170.197154'] == '0|1'].POS, otherSNPs[otherSNPs['IPF3170.197154'] == '0|1']['Sample Name'],
        mfc='orange', marker='s', mec='orangered', ms=8, linestyle='none', label='0|1')

ax1.plot(missenseQV[missenseQV['IPF3170.197154'] == '1|0'].POS, missenseQV[missenseQV['IPF3170.197154'] == '1|0']['Sample Name'],
        mfc='skyblue', marker='^', mec='blue', ms=8, linestyle='none', label='QV 1|0')

positions = np.arange(14500000, 14740000, 25000).tolist()
labels = ' '.join('{:,}'.format(x) for x in positions).split()
ax1.xaxis.set_major_locator(ticker.FixedLocator(positions))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.margins(0.01)
plt.grid(ls='--')
plt.subplots_adjust(left=0.07,
                    bottom=0.11,
                    right=0.98,
                    top=0.88,
                    wspace=0.2,
                    hspace=0.2)