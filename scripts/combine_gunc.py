import pandas as pd
import os.path

# Directing python to the input from snakemake
gunc_reports = snakemake.input["gunc_reports"]


dataframes = []
for file in gunc_reports:
    df = pd.read_csv(file+"/GUNC.progenomes_2.1.maxCSS_level.tsv", sep="\t")
    dataframes.append(df)

combined = pd.concat(dataframes)

combined = combined[
    [
        "genome",
        "n_genes_called",
        "n_genes_mapped",
        "n_contigs",
        "taxonomic_level",
        "proportion_genes_retained_in_major_clades",
        "genes_retained_index",
        "clade_separation_score",
        "contamination_portion",
        "n_effective_surplus_clades",
        "mean_hit_identity",
        "reference_representation_score",
        "pass.GUNC"
    ]
]

# combined.columns = [
#     "name",
#     "scientific_name",
#     "taxonomy_id",
#     "taxonomy_lvl",
#     "kraken_assigned_reads",
#     "added_reads",
#     "new_est_reads",
#     "fraction_total_reads",
# ]

combined.to_csv(snakemake.output[0], sep="\t", index=False)
#
