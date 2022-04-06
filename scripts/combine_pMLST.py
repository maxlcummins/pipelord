import pandas as pd
import os.path
import re

# Directing python to the input from snakemake
pmlst_summaries = snakemake.input["pMLST_summary"]

#Create an empty list
dataframes = []

for file in pmlst_summaries:
    #Read it in as df
    df = pd.read_csv(file, sep="\t", header=None, index_col=None)
    name = re.sub(".out/results_simple.txt","", file)
    name = re.sub(".*/", "", name)
    df['name'] = name
    scheme = re.sub(".*pMLST/","", file)
    scheme = re.sub("/.*", "", scheme)
    df['scheme'] = scheme
    dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Change the filename column to be pMLST
combined.rename(columns={0:'pMLST'}, inplace=True)

#Change the filename column to be pMLST
combined.rename(columns={'scheme':'pMLST_scheme'}, inplace=True)

#Change column order
combined = combined[[
    'name',
    'pMLST_scheme',
    'pMLST'
    ]
]

#Write our output file to text
combined.to_csv(snakemake.output[0], sep="\t", index=False)
#
