import pandas as pd
import os.path

# Directing python to the input from snakemake
mlst_summaries = snakemake.input["mlst_summary"]

#Create an empty list
dataframes = []

#For each of our mlst files
for file in mlst_summaries:
    #Read it in as df
    df = pd.read_csv(file, sep="\t", header=None, index_col=None)
    #Add the df to our list of dfs
    dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Extract the number of columns
numCols = combined.shape[1]

#Define the first three column names
mlst_colnames  =  [
    'name',
    'scheme',
    'ST'
]

#Create a list of column names for our alleles
allele_colnames = ["allele"+str(x) for x in range(1, numCols-2)]

#Combine our lists of column names
mlst_colnames = mlst_colnames + allele_colnames

#Assign our list of column names
combined.columns = mlst_colnames

#Remove the path related info from our name column
combined["name"] = combined["name"].str.replace(".*/", "", regex=True)
combined["name"] = combined["name"].str.replace(".fasta", "", regex=False)

#Fill empty columns with 0
combined = combined.fillna(0)

#Write our output file to text
combined.to_csv(snakemake.output[0], sep="\t", index=False)
#
