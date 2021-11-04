import pandas as pd
import os.path

# Directing python to the input from snakemake
bracken_reports = snakemake.input["bracken_reports"]


dataframes = []
for file in bracken_reports:
    df = pd.read_csv(file, sep="\t")
    df["source"] = os.path.basename(file)
    df["source"] = df["source"].str.replace(".bracken.txt", "", regex=False)
    dataframes.append(df)

combined = pd.concat(dataframes)

combined = combined[
    [
        "source",
        "name",
        "taxonomy_id",
        "taxonomy_lvl",
        "kraken_assigned_reads",
        "added_reads",
        "new_est_reads",
        "fraction_total_reads",
    ]
]

combined.columns = [
    "name",
    "scientific_name",
    "taxonomy_id",
    "taxonomy_lvl",
    "kraken_assigned_reads",
    "added_reads",
    "new_est_reads",
    "fraction_total_reads",
]

combined.to_csv(snakemake.output[0], sep="\t", index=False)
#
