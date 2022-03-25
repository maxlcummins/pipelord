import pandas as pd
import os.path

# Directing python to the inputs from snakemake
names_list = snakemake.params["names"]
gunc_report = snakemake.input["gunc_report"]
bracken_report = snakemake.input["bracken_report"]
assembly_stats_report = snakemake.input["assembly_stats_report"]
checkm_report = snakemake.input["checkm_report"]

#Create our df of sample names
names = df = pd.DataFrame (names_list, columns = ['name'])

#Read in our df
gunc = pd.read_csv(gunc_report, sep="\t", index_col=None)

#Prefix our columns
gunc = gunc.add_prefix('GUNC.')

#Rename some select columns
gunc.rename(columns={'GUNC.genome':'name', 'GUNC.pass.GUNC':'GUNC.pass'}, inplace=True)

gunc_simple = gunc[
    [
    'name',
    'GUNC.n_genes_called',
    'GUNC.n_genes_mapped',
    'GUNC.reference_representation_score',
    'GUNC.pass'
    ]
    ]

#Read in our df
bracken = pd.read_csv(bracken_report, sep="\t", index_col=None)

#Prefix our columns
bracken = bracken.add_prefix('BRACKEN.')

#Rename some select columns
bracken.rename(columns={'BRACKEN.name':'name'}, inplace=True)

#Group by name and filter for the species best match by proportion of reads mapped
top_species = bracken.groupby(['name'])['BRACKEN.fraction_total_reads'].transform(max) == bracken['BRACKEN.fraction_total_reads']
bracken_top_species = bracken[top_species]

#Select columns of interest
bracken_top_species = bracken_top_species[['name', 'BRACKEN.scientific_name', 'BRACKEN.fraction_total_reads']]

#Rename some select columns
bracken_top_species.rename(columns={'BRACKEN.scientific_name':'BRACKEN.species_best_match'}, inplace=True)

#Extract the second best match for species by proportion of reads mapped
not_max = bracken.groupby(['name'])['BRACKEN.fraction_total_reads'].transform(max) != bracken['BRACKEN.fraction_total_reads']
second_top_species = bracken[not_max]

second_top_species_idx = second_top_species.groupby(['name'])['BRACKEN.fraction_total_reads'].transform(max) == second_top_species['BRACKEN.fraction_total_reads']

#Extract the second most common match
bracken_second_top_species = second_top_species[second_top_species_idx]

#Select columns of interest
bracken_second_top_species = bracken_second_top_species[['name', 'BRACKEN.scientific_name', 'BRACKEN.fraction_total_reads']]

#Rename some select columns
bracken_second_top_species.rename(columns={'BRACKEN.scientific_name':'BRACKEN.species_second_best_match', 'BRACKEN.fraction_total_reads':'BRACKEN.fraction_total_reads_second_best_match'}, inplace=True)

#Combine our bracken dataframes
bracken_simple = pd.merge(bracken_top_species, bracken_second_top_species, left_on="name", right_on="name", how="outer")

print(bracken_simple)

#Read in our df
assembly_stats = pd.read_csv(assembly_stats_report, sep="\t", index_col=None)

#Prefix our columns
assembly_stats = assembly_stats.add_prefix('ASSEMBLY_STATS.')

#Rename some select columns
assembly_stats.rename(columns={'ASSEMBLY_STATS.name':'name'}, inplace=True)

#Select some columns of interest
assembly_stats_simple = assembly_stats[
    [
    'name',
    'ASSEMBLY_STATS.total_length',
    'ASSEMBLY_STATS.number',
    'ASSEMBLY_STATS.N50'
    ]
    ]

#Read in our df
checkm = pd.read_csv(checkm_report, sep="\t", index_col=None)

#Prefix our columns
checkm = checkm.add_prefix('CHECKM.')

#Rename some select columns
checkm.rename(columns={'CHECKM.Bin Id':'name'}, inplace=True)

#Select some columns of interest
checkm_simple = checkm[
    [
    'name',
    'CHECKM.Completeness',
    'CHECKM.Contamination',
    'CHECKM.Strain heterogeneity'
    ]
    ]

#Initialise our empty QC report df
summary_df = names

if path.exists("resources/test_data/test_output/QC_workflow/summaries/bracken_report.txt"):
    summary_df = pd.merge(summary_df, bracken_second_top_species, left_on="name", right_on="name", how="outer")
if path.exists("resources/test_data/test_output/QC_workflow/summaries/gunc_report.txt"):
    summary_df = pd.merge(summary_df, gunc_simple, left_on="name", right_on="name", how="outer")
if path.exists("resources/test_data/test_output/QC_workflow/summaries/checkm_qa.txt"):
    summary_df = pd.merge(summary_df, checkm_simple, left_on="name", right_on="name", how="outer")
if path.exists("resources/test_data/test_output/QC_workflow/summaries/assembly_stats.txt"):
    summary_df = pd.merge(summary_df, assembly_stats_simple, left_on="name", right_on="name", how="outer")

combined.to_csv(snakemake.output[0], sep="\t", index=False)
#