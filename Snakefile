import re
import os

configfile: "misc/masterconfig4.yaml"

outdir = config['outdir']
prefix = config['prefix']
gene_dbs = expand(config['gene_dbs'])
scheme = config['pmlst_scheme']

if path.exists("tools") == False:
    print('tools directory not located, creating tools directory...')
    os.system('mkdir tools')

if config['input_type'] == 'raw_reads':
    sample_ids, = glob_wildcards(config['raw_reads_path']+"/{sample}.R1.fastq.gz")
    include: "smks/read_cleaning.smk"
elif config['input_type'] == 'reads':
    sample_ids, = glob_wildcards(config['reads_path']+"/{sample}.R1.fastq.gz")
elif config['input_type'] == 'assemblies':
    sample_ids, = glob_wildcards(config['assembly_path']+"/{sample}.fasta")
else:
    print("\n")
    print('Please provide an input type of either raw_reads, reads or assemblies and ensure that you have provided the path to these files')
    print("\n")
    sys.exit()


#print(sample_ids,)
if len([i for i in sample_ids if '.' in i]) > 0:
    print("\n\n")
    print("Warning: One or more of your read sets has a full-stop/period (\'.\') in some of the sample names, shown below.")
    print('This is a problem because this character breaks many bioinformatic tools.')
    print('It is a character usually reserved for separating sample names from suffixes like \'.fasta\'.')
    print([i for i in sample_ids if '.' in i])
    print("\n\n")
    sys.exit()

#if any(r.match(sample_ids) for sample_ids in sample_ids):
 #   print(r.match(sample_ids) for sample_ids in sample_ids)

rule all:
    input:
        expand(config['outdir']+"/{prefix}/summaries/mlst.txt", prefix=prefix),
        expand(config['outdir']+"/{prefix}/summaries/Pointfinder.txt", prefix=prefix),
        expand(config['outdir']+"/{prefix}/abricate/{gene_db}/{sample}.tab", sample=sample_ids, gene_db=gene_dbs, prefix=prefix),
        expand(config['outdir']+"/{prefix}/summaries/abricate_hits.txt", prefix=prefix),
        expand(config['outdir']+"/{prefix}/summaries/pMLST.txt", prefix=config['prefix']),
        expand(config['outdir']+"/{prefix}/summaries/kraken2_full_summary.txt", prefix = config['prefix']),
        expand(config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta", sample=sample_ids, prefix=config['prefix']),
        # Summaries
        #expand(config['outdir']+"/{prefix}/summaries/fastp_summary.json", prefix=config['prefix']),
        expand(config['outdir']+"/{prefix}/summaries/assembly_stats.txt", prefix = config['prefix']),
        

#if config["general"]["seq_rep"] == "OTU" else [],

include: "smks/genome_assembly.smk"
include: "smks/species_id.smk"
include: "smks/genotype_abricate.smk"
include: "smks/point_mutations.smk"
include: "smks/strain_mlst.smk"
include: "smks/plasmid_mlst.smk"


#ruleorder: concatenate_abricate_hits > abricate_summary