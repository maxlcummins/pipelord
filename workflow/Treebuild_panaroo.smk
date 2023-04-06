import re
import subprocess
import os
from os import path
import git
import re

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

if config["genotype_modules"]['genome_annotater'] == "dfast":
    #For each of the subsets of samples we want to run
    for subset_prefix in subset_prefixes:
        #Report on the prefix name
        print("Subset with prefix \'"+subset_prefix+"\' contains the following genomes:")
        #Open the associated file of file names
        with open(subset_fofn[i]) as f:
            subset_filenames = ([line.strip("\n") for line in f if line.strip() != ''])
            #Print the list of genomes in the subset (from the file of file names)
            print(subset_filenames)
            #Path to gffs
            gff_path = [config['outdir']+"/"+config['prefix']+"/dfast/gffs/"+name+".gff" for name in subset_filenames]
            #Generate a list of files for our subsets which will be used by panaroo
            subsets_filenames.append(gff_path)
        #Increment our counter to help us dynamically access our files of file names
        i += 1
elif config["genotype_modules"]['genome_annotater'] == "prokka":
        #For each of the subsets of samples we want to run
    for subset_prefix in subset_prefixes:
        #Report on the prefix name
        print("Subset with prefix \'"+subset_prefix+"\' contains the following genomes:")
        #Open the associated file of file names
        with open(subset_fofn[i]) as f:
            subset_filenames = ([line.strip("\n") for line in f if line.strip() != ''])
            #Print the list of genomes in the subset (from the file of file names)
            print(subset_filenames)
            #Path to gffs
            gff_path = [config['outdir']+"/"+config['prefix']+"/prokka/gffs/"+name+".gff" for name in subset_filenames]
            #Generate a list of files for our subsets which will be used by panaroo
            subsets_filenames.append(gff_path)
        #Increment our counter to help us dynamically access our files of file names
        i += 1

print(subsets_filenames)

print(subset_prefixes)

rule all:
    input:
        expand(config['outdir']+"/{prefix}/subsets/{subset_prefix}_gff_locations.txt", prefix=prefix, subset_prefix=subset_prefixes),
        expand(config['outdir']+"/{prefix}/panaroo/{subset_prefix}/core_gene_alignment.aln", prefix=prefix, subset_prefix=subset_prefixes),
        expand(config['outdir']+"/{prefix}/panaroo_QC/{subset_prefix}", prefix=prefix, subset_prefix=subset_prefixes) if config["treebuild_modules"]["run_panaroo_QC"] else [],
        expand(config['outdir']+"/{prefix}/snp_sites/{subset_prefix}/core_gene_alignment_snp_sites.aln", prefix=prefix, subset_prefix=subset_prefixes),
        expand(config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_snp_sites.treefile", prefix=prefix, subset_prefix=subset_prefixes) if config["treebuild_modules"]["run_core_genome_treebuild_snp_sites"] else [],
        expand(config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_full.treefile", prefix=prefix, subset_prefix=subset_prefixes) if config["treebuild_modules"]["run_core_genome_treebuild_full"] else [],
        expand(config['outdir']+"/{prefix}/snp_dists/{subset_prefix}/CGA_snp_sites_pairwise_snps_counts.csv", prefix=prefix, subset_prefix=subset_prefixes) if config["treebuild_modules"]["run_panaroo"] else [],
        expand(config['outdir']+"/{prefix}/snp_dists/{subset_prefix}/CGA_full_pairwise_snps_counts.csv", prefix=prefix, subset_prefix=subset_prefixes) if config["treebuild_modules"]["run_panaroo"] else [],

if config["genotype_modules"]['genome_annotater'] == "dfast":
    rule make_fofns:
        input: config['subset_fofn'],
        output: expand(config['outdir']+"/{prefix}/subsets/{subset_prefix}_gff_locations.txt", prefix=prefix, subset_prefix=subset_prefixes)
        run:
            for f, o in zip(input,output):
                with open(f, 'r') as f_in:
                    with open(o, 'w') as f_out:
                        for line in f_in:
                            line_fix = re.sub('^', config['outdir']+"/"+config['prefix']+"/dfast/gffs/", line)
                            line_fix = re.sub('\n', '', line_fix)
                            line_fix = re.sub('$', ".gff\n", line_fix)
                            line_fix = re.sub(config['outdir']+"/"+config['prefix']+"/dfast/gffs/"+'.gff\n', "", line_fix)
                            print(line_fix)
                            f_out.write(line.replace(line, line_fix))
elif config["genotype_modules"]['genome_annotater'] == "prokka":
    rule make_fofns:
        input: config['subset_fofn'],
        output: expand(config['outdir']+"/{prefix}/subsets/{subset_prefix}_gff_locations.txt", prefix=prefix, subset_prefix=subset_prefixes)
        run:
            for f, o in zip(input,output):
                with open(f, 'r') as f_in:
                    with open(o, 'w') as f_out:
                        for line in f_in:
                            line_fix = re.sub('^', config['outdir']+"/"+config['prefix']+"/prokka/gffs/", line)
                            line_fix = re.sub('\n', '', line_fix)
                            line_fix = re.sub('$', ".gff\n", line_fix)
                            line_fix = re.sub(config['outdir']+"/"+config['prefix']+"/prokka/gffs/"+'.gff\n', "", line_fix)
                            print(line_fix)
                            f_out.write(line.replace(line, line_fix))
elif config["genotype_modules"]['genome_annotater'] == "bakta":
    rule make_fofns:
        input: config['subset_fofn'],
        output: expand(config['outdir']+"/{prefix}/subsets/{subset_prefix}_gff_locations.txt", prefix=prefix, subset_prefix=subset_prefixes)
        run:
            for f, o in zip(input,output):
                with open(f, 'r') as f_in:
                    with open(o, 'w') as f_out:
                        for line in f_in:
                            line_fix = re.sub('^', config['outdir']+"/"+config['prefix']+"/bakta/gffs/", line)
                            line_fix = re.sub('\n', '', line_fix)
                            line_fix = re.sub('$', ".gff\n", line_fix)
                            line_fix = re.sub(config['outdir']+"/"+config['prefix']+"/bakta/gffs/"+'.gff3\n', "", line_fix)
                            print(line_fix)
                            f_out.write(line.replace(line, line_fix))
else: print("Please set the genome annotater in the config file to either 'dfast', 'bakta' or 'prokka'")

rule panaroo_mash_db_dl:
    output:
        "resources/panaroo/refseq.genomes.k21s1000.msh"
    conda:
        "envs/panaroo.yaml"
    log:
        config['base_log_outdir']+"/"+config['prefix']+"/panaroo/QC_mash_db_dl.log"
    threads: 1
    shell:
        """
        wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh
        mv refseq.genomes.k21s1000.msh {output}
        """


rule panaroo_QC:
    input:
        gff_locs = config['outdir']+"/{prefix}/subsets/{subset_prefix}_gff_locations.txt",
        mash_db = "resources/panaroo/refseq.genomes.k21s1000.msh"
    output:
        directory(config['outdir']+"/{prefix}/panaroo_QC/{subset_prefix}")
    conda:
        "envs/panaroo.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/panaroo/{subset_prefix}_mashQC.log"
    threads: 12
    shell:
        """
        mkdir -p {output}
        cat {input.gff_locs} | xargs panaroo-qc -o {output} -t {threads} --graph_type all --ref_db {input.mash_db} -i
        """


rule panaroo:
    input:
        mashQC = config['outdir']+"/{prefix}/panaroo_QC/{subset_prefix}" if config["treebuild_modules"]["run_panaroo_QC"] else [],
        gff_locs = config['outdir']+"/{prefix}/subsets/{subset_prefix}_gff_locations.txt"
    output:
        outdir = directory(config['outdir']+"/{prefix}/panaroo/{subset_prefix}"),
        aln = config['outdir']+"/{prefix}/panaroo/{subset_prefix}/core_gene_alignment.aln"
    conda:
        "envs/panaroo.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/panaroo/{subset_prefix}.log"
    threads: 12
    params:
        out = config['outdir']+"/{prefix}/panaroo/{subset_prefix}",
        core_req = config['panaroo']['core_req'],
        aligner = config['panaroo']['aligner'],
        strictness = config['panaroo']['strictness']
    shell:
        "cat {input.gff_locs} | xargs panaroo -o {output.outdir} --clean-mode {params.strictness} -a core --aligner {params.aligner} --core_threshold {params.core_req} -t {threads} -i"

rule snp_sites:
    input:
        config['outdir']+"/{prefix}/panaroo/{subset_prefix}/core_gene_alignment.aln"
    output:
        config['outdir']+"/{prefix}/snp_sites/{subset_prefix}/core_gene_alignment_snp_sites.aln"
    conda:
        "envs/snp_sites.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/panaroo/logs/{subset_prefix}/snp_sites.log"
    shell:
        "snp-sites -c {input} > {output}"

rule iqtree:
    input:
        config['outdir']+"/{prefix}/panaroo/{subset_prefix}/core_gene_alignment.aln"
    output:
        config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_full.treefile"
    conda:
        "envs/iqtree.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/panaroo/logs/{subset_prefix}/iqtree.log"
    params:
        out = config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_full"
    threads: 12
    shell:
        """
        iqtree -s {input} -pre {params.out} -m GTR+F -ntmax {threads} -bb 1000
        """

rule iqtree_snp_sites:
    input:
        config['outdir']+"/{prefix}/snp_sites/{subset_prefix}/core_gene_alignment_snp_sites.aln"
    output:
        config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_snp_sites.treefile"
    conda:
        "envs/iqtree.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/panaroo/logs/{subset_prefix}/iqtree_snp_sites.log"
    params:
        out = config['outdir']+"/{prefix}/trees/{subset_prefix}/CGA_snp_sites"
    threads: 12
    shell:
        """
        iqtree -s {input} -pre {params.out} -m GTR+ASC -ntmax {threads} -bb 1000
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
            config['outdir']+"/{prefix}/panaroo/{subset_prefix}"
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


if config["treebuild_modules"]["run_panaroo_QC"] and config["treebuild_modules"]["run_core_genome_treebuild_snp_sites"] and config["treebuild_modules"]["run_core_genome_treebuild_full"]:
    ruleorder: make_fofns > panaroo_mash_db_dl > panaroo_QC > panaroo > snp_sites > iqtree_snp_sites > snp_dists_snp_sites > iqtree > snp_dists_full

if config["treebuild_modules"]["run_panaroo_QC"] and config["treebuild_modules"]["run_core_genome_treebuild_snp_sites"] == False and config["treebuild_modules"]["run_core_genome_treebuild_full"]:
    ruleorder: make_fofns > panaroo_mash_db_dl > panaroo_QC > panaroo > snp_sites > iqtree_snp_sites > snp_dists_snp_sites

if config["treebuild_modules"]["run_panaroo_QC"] and config["treebuild_modules"]["run_core_genome_treebuild_snp_sites"] and config["treebuild_modules"]["run_core_genome_treebuild_full"] == False:
    ruleorder: make_fofns > panaroo_mash_db_dl > panaroo_QC > panaroo > snp_sites > iqtree > snp_dists_full
