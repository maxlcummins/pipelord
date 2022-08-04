import pandas as pd
import os.path
import re

# Directing python to the input from snakemake
pmlst_summaries = snakemake.input["pMLST_summary"]

#Create an empty list
dataframes = []

#Initialise a for loop for processing our pMLST outputs
for file in pmlst_summaries:
    #Read it in as df
    df = pd.read_csv(file, sep="\t", header=None, on_bad_lines='skip')
    #Trim our filenames to use as a name column
    name = re.sub(".out/results.txt","", file)
    name = re.sub(".*/", "", name)
    #Assign our name column
    df['name'] = name
    #Reassign name columns. The input isnt actually a CSV so we label the original input as data and redundantly reassign our name column for the names
    df.columns = ['datacol', 'name']
    #Isolate pMLST scheme info
    df_profile = df[df.datacol.str.contains('pMLST profile')]
    #Assign colnames to our new df with pMLST scheme info
    df_profile.columns = ['pMLST_scheme', 'name']
    #Isolate pMLST info    
    df_pMLST = df[df.datacol.str.contains('Sequence Type')]
    #Assign colnames to our new df with pMLST info
    df_pMLST.columns = ['Sequence_type', 'name']
    #Merge our dataframes with this new info
    df2 = pd.merge(df_profile, df_pMLST)
    #Check if there are any lines which indicate that no MLST loci were found
    no_loci = df.datacol.str.contains('No MLST loci was found').any()
    #Check if there are any lines which indicate that either a novel ST was found or no ST was found
    #pMLST tool doesn't discriminate these two particularly well on this line.
    unknown_pMLST = df.datacol.str.contains('Unknown').any()
    #If the file indicates no MLST loci were found AND that the genome has a novel or NA ST then assign it as 'None'
    if no_loci > 0 and unknown_pMLST > 0:
        df2.Sequence_type = 'None'
    #If the file indicates that some MLST loci were found but the tool still couldnt assign a pMLST then assign it as 'Novel'
    elif no_loci == 0 and unknown_pMLST > 0:
        df2.Sequence_type = 'Novel'
    #Pull out lines which indicate the closest STs
    df_nearest_pMLST = df[df.datacol.str.contains('Nearest ST')]
    #Assign column names to our new dataframe with nearest STs
    df_nearest_pMLST.columns = ['Nearest_ST', 'name']
    #Trim unwanted chars preceeding nearest ST data
    df_nearest_pMLST.Nearest_ST = df_nearest_pMLST['Nearest_ST'].str.replace("Nearest ST: ","", regex=False)
    #Trim unwanted chars preceeding nearest ST data. Same again but for when multiple matches of equal weight are identified
    df_nearest_pMLST.Nearest_ST = df_nearest_pMLST['Nearest_ST'].str.replace("Nearest STs: ","", regex=False)
    #Merge the dataframe with nearest ST data with our info on scheme and ST
    df2 = pd.merge(df2, df_nearest_pMLST, on='name', how='outer')
    #Trim unwanted chars preceeding scheme info
    df2.pMLST_scheme = df2.pMLST_scheme.str.replace("pMLST profile: ","", regex=False)
    #Trim unwanted chars - after a space is just the word pmlst
    df2.pMLST_scheme = df2.pMLST_scheme.str.replace(" .*","", regex=True)
    #Trim unwanted chars in IncA/C scheme - slashes are naughty characters
    df2.pMLST_scheme = df2.pMLST_scheme.str.replace("/","", regex=False)
    #Convert schemes lowercase
    df2.pMLST_scheme = df2.pMLST_scheme.str.lower()
    #Trim unwanted chars preceeding ST
    df2.Sequence_type = df2.Sequence_type.str.replace("Sequence Type: ","", regex=False)
    #Trim unwanted chars before and after our IncF RSTs
    df2.Sequence_type = df2.Sequence_type.str.replace("[","", regex=False)
    df2.Sequence_type = df2.Sequence_type.str.replace("]","", regex=False)
    #Append our dataframes
    dataframes.append(df2)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Change column order
combined = combined[[
    'name',
    'pMLST_scheme',
    'Sequence_type',
    'Nearest_ST'
    ]
]
#Write our output file to text
combined.to_csv(snakemake.output[0], sep="\t", index=False)
#
