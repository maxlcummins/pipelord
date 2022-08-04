import pandas as pd
import os.path
import re

# Directing python to the input from snakemake
pmlst_summaries = snakemake.input['pMLST_summary']

#Create an empty list
dataframes = []

#Initialise a for loop for processing our pMLST outputs
for file in pmlst_summaries:
    #Read it in as df
    df = pd.read_csv(file, sep="\t", header=None, on_bad_lines='skip')
    #Trim our filenames to use as a name column
    name = re.sub(".out/results.txt","", file)
    name = re.sub(".*/", "", name)
    #Assign a column name to our input which isnt actually a CSV and therefore has a heading of '0'
    df.columns = ['datacol']
    #Identify pMLST data pertaining to IncF plasmids
    IncF_results = df['datacol'].str.contains('IncF RST').any()
    #Only process IncF data.
    #Wasteful loop iterations but snakemake was being painful so I decided on feeding all pMLST input rather than just IncF
    if IncF_results:
        #Pull out the columns with F replicons
        df = df[df.datacol.str.contains('^FI')]
        #Convert multiple spaces to commas to properly delimit our file
        df['datacol'] = df['datacol'].str.replace('  +', ',', regex=True)
        #Split our string by comma delimiter into different columns
        df = df['datacol'].str.split(",", expand=True)
        #Assign our name column
        df['name'] = name
        #Append our dataframes
        dataframes.append(df)

#Concatenate our dfs into a single one
combined = pd.concat(dataframes)

#Change column order
combined.columns = [
    'Locus',
    'Identity',
    'Coverage',
    'Alignment Length',
    'Allele Length',
    'Gaps',
    'Allele',
    'name'
    ]

#Change column order
combined = combined[[
    'name',
    'Locus',
    'Identity',
    'Coverage',
    'Alignment Length',
    'Allele Length',
    'Gaps',
    'Allele'
    ]
]

#Simplify when no allele is found with a dash
combined['Allele'] = combined['Allele'].str.replace('No hit found', '-', regex=False)

#Create a column for a simplified allele call
combined['Simple Allele'] = combined['Allele']
combined['Simple Allele'] = combined['Simple Allele'].str.replace('FII[S|Y|K]', '', regex=True)
combined['Simple Allele'] = combined['Simple Allele'].str.replace('FI', '', regex=True)
combined['Simple Allele'] = combined['Simple Allele'].str.replace('I', 'F', regex=True)
combined['Simple Allele'] = combined['Simple Allele'].str.replace('_', '', regex=True)

#Separate our alleles to combine them in the correct order for IncF RST
#Sort by allele name. If they are instead sorted in order of discovery (via contig number) then we will see inconsitencies between samples
FICs = combined[combined['Locus'] == 'FIC'].sort_values(by='Simple Allele', ascending=False)
FIIs = combined[combined['Locus'] == 'FII'].sort_values(by='Simple Allele', ascending=False)
FIKs = combined[combined['Locus'] == 'FIIK'].sort_values(by='Simple Allele', ascending=False)
FISs = combined[combined['Locus'] == 'FIIS'].sort_values(by='Simple Allele', ascending=False)
FIYs = combined[combined['Locus'] == 'FIIY'].sort_values(by='Simple Allele', ascending=False)

FIA_type_df = combined[combined['Locus'] == 'FIA'].sort_values(by='Simple Allele', ascending=False)

FIB_type_df = combined[combined['Locus'] == 'FIB'].sort_values(by='Simple Allele', ascending=False)

#Combine our dataframes again
F_type_df = pd.concat([FIIs, FICs, FIKs, FISs, FIYs])
F_type_df = F_type_df[F_type_df['Simple Allele'] != '-']

#Group by name and collapse our FIIs, FICs, FIKs, FISs and FIYs with a slash separator
F_type = F_type_df.groupby('name').apply(lambda x: '/'.join(x['Simple Allele'])).reset_index()

#Assign our column names
F_type.columns = ['name', 'F type']

#Group by name and collapse our FIBs with a slash separator
FIA_type = FIA_type_df.groupby('name').apply(lambda x: '/'.join(x['Simple Allele'])).reset_index()

#Assign our column names
FIA_type.columns = ['name', 'FIA type']

#Group by name and collapse our FIBs with a slash separator
FIB_type = FIB_type_df.groupby('name').apply(lambda x: '/'.join(x['Simple Allele'])).reset_index()

#Assign our column names
FIB_type.columns = ['name', 'FIB type']

#Join our columns
IncF_RSTs = pd.merge(F_type, FIA_type, on='name',how='outer')
IncF_RSTs = pd.merge(IncF_RSTs, FIB_type, on='name',how='outer')

#Replace NaNs (for missing alleles) with a dash
IncF_RSTs = IncF_RSTs.fillna('-')

#Replace the cell contents for F types to make them better reflect null values 
IncF_RSTs['F type'] = IncF_RSTs['F type'].str.replace('-', 'F-', regex=True)
IncF_RSTs['FIA type'] = IncF_RSTs['FIA type'].str.replace('-', 'A-', regex=True)
IncF_RSTs['FIB type'] = IncF_RSTs['FIB type'].str.replace('-', 'B-', regex=True)

#Generate an IncF RST column by combining the alleles
IncF_RSTs['IncF RST'] = IncF_RSTs['F type']+':'+IncF_RSTs['FIA type']+':'+IncF_RSTs['FIB type']

#Make a simple IncF RST dataframe to bind to our full IncF RST dataframe
IncF_RSTs_simple = IncF_RSTs[['name', 'IncF RST']]

#Join our columns
combined = pd.merge(combined, IncF_RSTs_simple, on='name',how='outer')

#Write our simple file to text
IncF_RSTs.to_csv(snakemake.output[0], sep="\t", index=False)

#Write our detailed file to text
combined.to_csv(snakemake.output[1], sep="\t", index=False)

