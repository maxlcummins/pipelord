import pandas as pd
import os.path

# Directing python to the input from snakemake
assembly_stats_summaries = snakemake.input["pointfinder_summary"]
#Create an empty list
dataframes = []

#For each of our mlst files
for file in assembly_stats_summaries:
    #Read it in as df
    df = pd.read_csv(file, sep="\t", index_col=None)
    df["name"] = df["name"].str.replace("_blastn_results.tsv", "", regex=False)
    df["name"] = df["name"].str.replace(".*/", "", regex=True)
    #Add the df to our list of dfs
    dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Write our output file to text
combined.to_csv(snakemake.output[0], sep="\t", index=False)
#
