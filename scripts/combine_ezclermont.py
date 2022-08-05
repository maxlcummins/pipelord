import pandas as pd
import os.path
import re

# Directing python to the input from snakemake
ezclermont_summaries = snakemake.input['ezclermont_summary']

print(ezclermont_summaries)

#Create an empty list
dataframes = []

#Initialise a for loop for processing our pMLST outputs
for file in ezclermont_summaries:
    #Read it in as df
    df = pd.read_csv(file, sep="\t", header=None, on_bad_lines='skip')
    #Assign a column name to our input which isnt actually a CSV and therefore has a heading of '0'
    df.columns = ['datacol']
    name = df[df.datacol.str.contains('Reading in sequence')]  
    #Trim our filenames to use as a name column
    name.datacol = name.datacol.str.replace('Reading in sequence(s) from', '', regex=False)
    name.datacol = name.datacol.str.replace('.fasta', '', regex=False)
    #Create an empty dataframe
    df2 = pd.DataFrame()
    #Create our name column
    df2['name'] = name
    #Save the result for trpBA
    df2['trpBA'] = df[df.datacol.str.contains('trpBA_control')].iloc[0,0]
    print(df2['trpBA'])
    #Save the result for virA
    df2['virA'] = df[df.datacol.str.contains('virA:')].iloc[0,0]
    #Save the result for TspE4
    df2['TspE4'] = df[df.datacol.str.contains('TspE4:')].iloc[0,0]
    #Save the result for arpA:
    df2['arpA'] = df[df.datacol.str.contains('arpA:')].iloc[0,0]
    #Save the result for chu:
    df2['chu'] = df[df.datacol.str.contains('chu:')].iloc[0,0]
    #Save the result for yjaA:
    df2['yjaA'] = df[df.datacol.str.contains('yjaA:')].iloc[0,0]
    #Save the result for arpA:
    df2['arpA'] = df[df.datacol.str.contains('arpA:')].iloc[0,0]
    #Save the result for the Clermont_type:
    df2['Clermont_type'] = df[df.datacol.str.contains('Clermont type:')].iloc[0,0]
    #Append our dataframes
    dataframes.append(df2)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Select our columns we wish to change
colnames = combined.select_dtypes(['object']).columns

#Replace unwanted characters in our dataframe
combined[colnames] = combined[colnames].replace('.*: ', '', regex=True)

#Replace dashes and pluses with 1s and 0s
combined[colnames] = combined[colnames].replace('+', 1, regex=False)
combined[colnames] = combined[colnames].replace('-', 0, regex=False)

#Write our simple file to text
combined.to_csv(snakemake.output[0], sep="\t", index=False)

