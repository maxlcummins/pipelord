import pandas as pd
import os.path

# Directing python to the input from snakemake
fimtyper_summaries = snakemake.input["fimtyper_summary"]

#Read in out dataframe and add a name column
dataframes = []
for file in fimtyper_summaries:
    df = pd.read_csv(file, sep="\t")
    if 'Contig'in df.columns:
        ## If detect No FimH type found skip sample
        ## If detect "Please contact curator, if new fimH-type" skip line
        df["name"] = os.path.basename(file)
        df["name"] = df["name"].str.replace("_named.txt", "", regex=False)
        dataframes.append(df)
    else:
        print("No fimH type found in "+str(os.path.basename(file)))

#Combine our fimtyper files
combined = pd.concat(dataframes)

#Name our columns
combined = combined[
    [
        "name",
        "Fimtype",
        "Identity",
        "Query/HSP",
        "Contig",
        "Position in contig",
        "Accession no."
    ]
]

#Reorder our columns (probably redundant now)
combined = combined[ ['name'] + [ col for col in combined.columns if col != 'name' ] ]

#Remove accession number column (usually empty)
combined = combined.drop('Accession no.', axis=1)

#Add a new column for query length and HSP (High-scoring segment pair i.e. match) length
combined["Query length"] = combined["Query/HSP"].str.replace("/.*", "", regex=True)
combined["HSP length"] = combined["Query/HSP"].str.replace(".*/", "", regex=True)

#Add an asterisk to novel alleles which have a non 100% identity or a length other than 489
combined.loc[combined['Identity'] != 100.0, 'Fimtype'] = combined.loc[combined['Identity'] != 100.0, 'Fimtype']+"*"
combined.loc[combined['Query length'] != combined['HSP length'], 'Fimtype'] = combined.loc[combined['Query length'] != combined['HSP length'], 'Fimtype']+"*"

#Remove instances of double asterisks
combined["Fimtype"] = combined["Fimtype"].str.replace("**", "*", regex=False)

#Define pattern of lines we want to remove (Messages from devs to report novel alleles)
false_line_pattern = 'Please'

#Identify rows which are lines we want to remove
pattern_filter = combined['Fimtype'].str.contains(false_line_pattern)
 
#Remove unwanted rows
combined = combined[~pattern_filter]

#Write to file
combined.to_csv(snakemake.output[0], sep="\t", index=False)

