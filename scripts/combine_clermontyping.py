import pandas as pd
import os.path
import re

# Directing python to the input from snakemake
clermontyping_summaries = snakemake.input['clermontyping']

print(clermontyping_summaries)

#Create an empty list
dataframes = []

#Initialise a for loop for processing our pMLST outputs
for file in clermontyping_summaries:
    #Read it in as df
    df = pd.read_csv(file, sep="\t", header=None)
    #Assign a column name to our input which isnt actually a CSV and therefore has a heading of '0'
    df.columns = ['name', 'genes_detected', 'quadraplex_resuts', 'alleles_for_C_or_G', 'clermontyping_phylogroup', 'mashfile']
    #Append our dataframes
    dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Replace unwanted characters in our dataframe
combined['name'] = combined['name'].replace('\.fasta', '', regex=True)

#Write our simple file to text
combined.to_csv(snakemake.output[0], sep="\t", index=False)

