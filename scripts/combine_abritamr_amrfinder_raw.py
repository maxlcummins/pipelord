import pandas as pd
import os.path
import re

# Directing python to the input from snakemake
amrfinder = snakemake.input["amrfinder"]

#Create an empty list
dataframes = []

for file in amrfinder:
    df = pd.read_csv(file, sep="\t", index_col=None)
    df["name"] = file
    df["name"] = df["name"].str.replace("/amrfinder.out", "", regex=False)
    df["name"] = df["name"].str.replace(".*\/", "", regex=True)
    dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Move the name column to the start
combined = combined[['name'] + [ col for col in combined.columns if col != 'name' ]]

#Write our output file to text
combined.to_csv(snakemake.output["abritamr_resistance"], sep="\t", index=False)
