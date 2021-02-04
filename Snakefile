configfile: "misc/masterconfig.yaml"

outdir = config['outdir']
prefix = config['prefix']
gene_dbs = expand(config['gene_dbs'])
sample_ids, = glob_wildcards(config['raw_reads_path']+"/{sample}.R1.fastq.gz")

rule all:
    input:
        expand(config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta", sample=sample_ids, prefix=config['prefix']),
        expand(config['outdir']+"/{prefix}/abricate/{gene_db}/{sample}.tab", sample=sample_ids, gene_db=gene_dbs, prefix=prefix),
        expand(config['outdir']+"/{prefix}/summaries/abricate_hits.txt", prefix=prefix),
        # Summaries
        expand(config['outdir']+"/{prefix}/summaries/fastp_summary.json", prefix=config['prefix']),
        expand(config['outdir']+"/{prefix}/summaries/assembly_stats.txt", prefix = config['prefix']),
        expand(config['outdir']+"/{prefix}/summaries/kraken2_summary.txt", prefix = config['prefix'])

include: "smks/read_cleaning.smk"
include: "smks/genome_assembly.smk"
include: "smks/species_id.smk"
include: "smks/genotype_abricate.smk"
