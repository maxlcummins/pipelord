import pandas as pd
import os.path

# Directing python to the input from snakemake
spifinder_summaries = snakemake.input["fimtyper_summary"]


dataframes = []
for file in spifinder_summaries:
    df = pd.read_csv(file, sep="\t")
    if 'Contig'in df.columns:
        ## If detect No FimH type found skip sample
        ## If detect "Please contact curator, if new fimH-type" skip line
        df["name"] = os.path.basename(file)
        df["name"] = df["name"].str.replace("_named.txt", "", regex=False)
        dataframes.append(df)
    else:
        print("No fimH type found in "+str(os.path.basename(file)))

combined = pd.concat(dataframes)

combined = combined[
    [
        "name",
        "Fimtype",
        "Identity",
        "Query/HSP",
        "Contig",
        "Position in contig",
        "Accession no."
    ]
]

combined = combined[ ['name'] + [ col for col in combined.columns if col != 'name' ] ]

combined = combined.drop('Accession no.', axis=1)

combined.to_csv(snakemake.output[0], sep="\t", index=False)
#
