import pandas as pd
import os.path
import re

#### Process our resistance data ####

# Directing python to the input from snakemake
abritamr_resistance = snakemake.input["abritamr_resistance"]

#Create an empty list
dataframes = []

for file in abritamr_resistance:
    df = pd.read_csv(file, sep="\t", index_col=None)
    dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Change the filename column to be name
combined.rename(columns={'Isolate':'name'}, inplace=True)

#Trim characters before the samplename remaining from the path
combined['name'].replace({r".*\/" : ""}, inplace = True, regex = True)

#Write our output file to text
combined.to_csv(snakemake.output["abritamr_resistance"], sep="\t", index=False)


#### Process our virulence data ####

# Directing python to the input from snakemake
abritamr_virulence = snakemake.input["abritamr_virulence"]

#Create an empty list
dataframes = []

for file in abritamr_virulence:
    df = pd.read_csv(file, sep="\t", index_col=None)
    dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Change the filename column to be name
combined.rename(columns={'Isolate':'name'}, inplace=True)

#Trim characters before the samplename remaining from the path
combined['name'].replace({r".*\/" : ""}, inplace = True, regex = True)

#Write our output file to text
combined.to_csv(snakemake.output["abritamr_virulence"], sep="\t", index=False)

#### Process our partial data ####

# Directing python to the input from snakemake
abritamr_partials = snakemake.input["abritamr_partials"]

#Create an empty list
dataframes = []

for file in abritamr_partials:
    df = pd.read_csv(file, sep="\t", index_col=None)
    dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Change the filename column to be name
combined.rename(columns={'Isolate':'name'}, inplace=True)

#Trim characters before the samplename remaining from the path
combined['name'].replace({r".*\/" : ""}, inplace = True, regex = True)

#Write our output file to text
combined.to_csv(snakemake.output["abritamr_partials"], sep="\t", index=False)