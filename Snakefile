configfile: "misc/masterconfig.yaml"

outdir = config['outdir']
prefix = config['prefix']
gene_dbs = expand(config['gene_dbs'])
scheme = config['pmlst_scheme']
sample_ids, = glob_wildcards(config['raw_reads_path']+"/{sample}.R1.fastq.gz")

rule all:
    input:
        expand(config['outdir']+"/{prefix}/summaries/mlst.txt", prefix=prefix),
        #expand(config['outdir']+"/{prefix}/summaries/{gene_db}_hits.txt", gene_db=gene_dbs, prefix=prefix),
        expand(config['outdir']+"/{prefix}/summaries/Pointfinder.txt", prefix=prefix),
        expand(config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta", sample=sample_ids, prefix=config['prefix']),
        expand(config['outdir']+"/{prefix}/abricate/{gene_db}/{sample}.tab", sample=sample_ids, gene_db=gene_dbs, prefix=prefix),
        expand(config['outdir']+"/{prefix}/summaries/abricate_hits.txt", prefix=prefix),
        expand(config['outdir']+"/{prefix}/summaries/pMLST.txt", prefix=config['prefix']),
        # Summaries
        expand(config['outdir']+"/{prefix}/summaries/fastp_summary.json", prefix=config['prefix']),
        expand(config['outdir']+"/{prefix}/summaries/assembly_stats.txt", prefix = config['prefix']),
        expand(config['outdir']+"/{prefix}/summaries/kraken2_full_summary.txt", prefix = config['prefix'])

#if config["general"]["seq_rep"] == "OTU" else [],

include: "smks/read_cleaning.smk"
include: "smks/genome_assembly.smk"
include: "smks/species_id.smk"
include: "smks/genotype_abricate.smk"
include: "smks/point_mutations.smk"
include: "smks/strain_mlst.smk"
include: "smks/plasmid_mlst.smk"


#ruleorder: concatenate_abricate_hits > abricate_summary