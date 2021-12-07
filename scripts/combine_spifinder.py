import pandas as pd
import os.path

# Directing python to the input from snakemake
spifinder_summaries = snakemake.input["spifinder_summary"]


dataframes = []
for file in spifinder_summaries:
    df = pd.read_csv(file, sep="\t")
    df["source"] = os.path.basename(file)
    df["source"] = df["source"].str.replace(".txt", "", regex=False)
    dataframes.append(df)

combined = pd.concat(dataframes)

combined = combined[
    [
        "Database",
        "SPI",
        "Identity",
        "Query / Template length",
        "Contig",
        "Position in contig",
        "Organism",
        "Insertion Site",
        "Category Function",
        "Note",
        "Accession number",
        "source"
    ]
]

combined.columns = [
    "Database",
    "SPI",
    "Identity",
    "Query / Template length",
    "Contig",
    "Position in contig",
    "Organism",
    "Insertion Site",
    "Category Function",
    "Note",
    "Accession number",
    "name"
]

combined = combined[ ['name'] + [ col for col in combined.columns if col != 'name' ] ]

combined.to_csv(snakemake.output[0], sep="\t", index=False)
#
