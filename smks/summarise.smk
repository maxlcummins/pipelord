import re
import subprocess
import os
from os import path
import git

rule summarise_all:
    input:
        mlst = expand(config['outdir']+"/{prefix}/summaries/mlst.txt", prefix=prefix),
        amr_snps = expand(config['outdir']+"/{prefix}/summaries/Pointfinder.txt", prefix=prefix),
        genotype = expand(config['outdir']+"/{prefix}/summaries/abricate_hits.txt", prefix=prefix),
        pmlst = expand(config['outdir']+"/{prefix}/summaries/pMLST.txt", prefix=prefix),
        #species_id = expand(config['outdir']+"/{prefix}/summaries/kraken2_full_summary.txt", prefix = prefix),
        #assembly_stats = expand(config['outdir']+"/{prefix}/summaries/assembly_stats.txt", prefix = prefix),
    output:
        expand(config['outdir']+"/{prefix}/summaries/{prefix}_simple_summary_N"+str(config['abricateR_identity'])+"L"+str(config['abricateR_length'])+".csv", prefix=prefix)
    conda:
        "config/R.yaml"
    params:
        prefix = config['prefix'],
        identity = config['abricateR_identity'],
        length = config['abricateR_length'],
        outdir = expand(config['outdir']+"/{prefix}/summaries/", prefix=prefix)

    shell:
        """
        #Rscript misc/abricateR.R --output {params.prefix} --output_directory {params.outdir} --abricate_in {input.genotype} --pointfinder_data {input.amr_snps} --pMLST_data {input.pmlst} --identity {params.identity} --length {params.length}
        Rscript misc/abricateR.R {input.genotype} {params.outdir} {params.prefix} {params.identity} {params.length} {input.amr_snps} {input.pmlst}
        """

#ruleorder: concatenate_abricate_hits > abricate_summary

#script misc/abricateR.R --abricate_in test_data/output/test/summaries/abricate_hits.txt --output test --output_directory test_data/output/test/summaries  --pointfinder_data test_data/output/test/summaries/Pointfinder.txt --pMLST_data test_data/output/test/summaries/pMLST.txt --identity 90 --length 90

#source("misc/abricateR.R")
#abricateR(abricate_in="test_data/output/test/summaries/abricate_hits.txt",output="test",output_directory="test_data/output/test/summaries",pointfinder_data="test_data/output/test/summaries/Pointfinder.txt",pMLST_data="test_data/output/test/summaries/pMLST.txt")

Rscript misc/abricateR.R "test" "test_data/output/test/summaries" "test_data/output/test/summaries/abricate_hits.txt" 90 90 "test_data/output/test/summaries/Pointfinder.txt" "test_data/output/test/summaries/pMLST.txt"