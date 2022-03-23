import pandas as pd
import os.path
import re

# Directing python to the input from snakemake
abricate_summaries = snakemake.input["abricate_summary"]

#Create an empty list
dataframes = []

for file in abricate_summaries:
    df = pd.read_csv(file, sep="\t", index_col=None)
    dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Change the filename column to be name
combined.rename(columns={'#FILE':'name'}, inplace=True)

#Write our output file to text
combined.to_csv(snakemake.output[0], sep="\t", index=False)
#
