import pandas as pd
import os.path

# Directing python to the input from snakemake
ectyper_summaries = snakemake.input["ectyper_summary"]

#Create an empty list
dataframes = []

#For each of our ectyper files
for file in ectyper_summaries:
    print(file)
    #Read it in as df
    df = pd.read_csv(file, sep="\t", header=0, index_col=None)
    #Add the df to our list of dfs
    dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Write our output file to text
combined.to_csv(snakemake.output[0], sep="\t", index=False)
#
