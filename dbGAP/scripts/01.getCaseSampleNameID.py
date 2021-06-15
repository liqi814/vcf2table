import pandas as pd

KIF15_variants = pd.read_csv("KIF15_variants/2021-06-07_RareEnsemble_KIF15_filter_variants_nlcluster.csv")
KIF15_variants[KIF15_variants.cohort == 'IPF'][['Indiv','VariantID']].to_csv("data/CaseSampleName_ID_dict_hg38.txt", sep='\t', header=False, index=False)