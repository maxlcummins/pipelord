import re
import subprocess
import os
from os import path
import git

configfile: "config/config_template.yaml"

#Define some parameters from our config file
outdir = config['outdir']
prefix = config['prefix']

email = config['email']

#Access our files of file names for our subsets
subset_fofn = config["subset_fofn"]

#Access our names for our subset groups
subset_prefixes = config["subset_name"]

#Initialise an empty list for later
subsets_filenames = []

#Set an incrementer for our for loops
i = 0

#For each of the subsets of samples we want to run
for subset_prefix in subset_prefixes:
    #Report on the prefix name
    print("Subset with prefix \'"+subset_prefix+"\' contains the following genomes:")
    #Open the associated file of file names
    with open(subset_fofn[i]) as f:
        subset_filenames = ([line.strip("\n") for line in f if line.strip() != ''])
        #Print the list of genomes in the subset (from the file of file names)
        print(subset_filenames)
        #Generate a list of files for our subsets which will be used by roary
        subsets_filenames.append(expand(config['outdir']+"/{prefix}/subsets/{subset_prefix}/gffs/{subset}.gff", prefix=prefix, subset_prefix=subset_prefix, subset=subset_filenames))
    #Increment our counter to help us dynamically access our files of file names
    i += 1

rule all:
    input:
        expand(config['outdir']+"/{prefix}/roary/{subset_prefix}", prefix=prefix, subset_prefix=subset_prefixes),
        expand(config['outdir']+"/{prefix}/snp_sites/{subset_prefix}/core_gene_alignment_snp_sites.aln", prefix=prefix, subset_prefix=subset_prefixes),
        expand(config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_snp_sites.treefile", prefix=prefix, subset_prefix=subset_prefixes) if config["treebuild_modules"]["run_core_genome_treebuild_snp_sites"] else [],
        expand(config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_full.treefile", prefix=prefix, subset_prefix=subset_prefixes) if config["treebuild_modules"]["run_core_genome_treebuild_full"] else [],
        expand(config['outdir']+"/{prefix}/snp_dists/{subset_prefix}/CGA_snp_sites_pairwise_snps_counts.csv", prefix=prefix, subset_prefix=subset_prefixes) if config["treebuild_modules"]["run_roary"] else [],
        expand(config['outdir']+"/{prefix}/snp_dists/{subset_prefix}/CGA_full_pairwise_snps_counts.csv", prefix=prefix, subset_prefix=subset_prefixes) if config["treebuild_modules"]["run_roary"] else [],


rule roary_subset:
    input:
        config['outdir']+"/{prefix}/dfast/gffs/{subset}.gff"
    output:
        config['outdir']+"/{prefix}/subsets/{subset_prefix}/gffs/{subset}.gff"
    shell:
        """
        cp -n {input} {output}
        """

rule roary:
    input:
        expand(config['outdir']+"/{prefix}/subsets/{subset_prefix}/gffs/{subset}.gff", prefix=prefix, subset_prefix=subset_prefix, subset=subset_filenames)
    output:
        directory(config['outdir']+"/{prefix}/roary/{subset_prefix}")
    conda:
        "envs/roary.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/roary/{subset_prefix}.log"
    threads:
        12
    params:
        out = config['outdir']+"/{prefix}/roary/{subset_prefix}",
        core_req = config['roary_core_req'],
        kraken_db = config['krakendb']
    shell:
        "roary -e -v -n -r -cd {params.core_req} -qc -k {params.kraken_db} -f {params.out} {input}  2> {log}"

rule snp_sites:
    input:
        config['outdir']+"/{prefix}/roary/{subset_prefix}"
    output:
        config['outdir']+"/{prefix}/snp_sites/{subset_prefix}/core_gene_alignment_snp_sites.aln"
    conda:
        "envs/snp_sites.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/roary/logs/{subset_prefix}/snp_sites.log"
    shell:
        "snp-sites -c {input}/core_gene_alignment.aln > {output}"

rule iqtree:
    input:
        config['outdir']+"/{prefix}/roary/{subset_prefix}"
    output:
        config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_full.treefile"
    conda:
        "envs/iqtree.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/roary/logs/{subset_prefix}/iqtree.log"
    params:
        out = config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_full"
    shell:
        """
        iqtree -s {input}/core_gene_alignment.aln -pre {params.out} -m GTR+F -bb 1000
        """

rule iqtree_snp_sites:
    input:
        config['outdir']+"/{prefix}/snp_sites/{subset_prefix}/core_gene_alignment_snp_sites.aln"
    output:
        config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_snp_sites.treefile"
    conda:
        "envs/iqtree.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/roary/logs/{subset_prefix}/iqtree_snp_sites.log"
    params:
        out = config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_snp_sites"
    shell:
        """
        iqtree -s {input} -pre {params.out} -m GTR+ASC -bb 1000
        """

rule snp_dists_snp_sites:
        input:
                config['outdir']+"/{prefix}/snp_sites/{subset_prefix}/core_gene_alignment_snp_sites.aln"
        output:
                config['outdir']+"/{prefix}/snp_dists/{subset_prefix}/CGA_snp_sites_pairwise_snps_counts.csv"
        conda:
                "envs/snp_dists.yaml"
        shell:
                "snp-dists -c {input} > {output}"



rule snp_dists_full:
        input:
                config['outdir']+"/{prefix}/roary/{subset_prefix}"
        output:
                config['outdir']+"/{prefix}/snp_dists/{subset_prefix}/CGA_full_pairwise_snps_counts.csv"
        conda:
                "envs/snp_dists.yaml"
        shell:
                "snp-dists -c {input}/core_gene_alignment.aln > {output}"


onerror:
    print("An error occurred")
    #Email user
    if config['email'] != '':
        shell("echo Your Snakemake job with prefix \'{prefix}\' has failed. | mail -s 'Snakemake job has failed' {email}")


onsuccess:
#Announce job finish
    print("Job finished")
    #Email user
    if config['email'] != '':
        shell("echo Your Snakemake job with prefix \'{prefix}\' has finished. It has been written to \'{outdir}/{prefix}/summaries/\' | mail -s 'Snakemake job has finished' {email}")


ruleorder: roary > snp_sites > iqtree > iqtree_snp_sites > snp_dists_snp_sites > snp_dists_full