import re
import subprocess
import os

prefix = config["prefix"]
outdir = config["outdir"]
maxthreads = snakemake.utils.available_cpu_count()

rule gunc_run:
	input:
		assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
	output:
		guncoutdir = directory(config["outdir"]+"/{prefix}/QC_workflow/gunc/{sample}.out"),
		guncout = config["outdir"]+"/{prefix}/QC_workflow/gunc/{sample}.out/GUNC.progenomes_2.1.maxCSS_level.tsv",
		tempdir = directory(temp(config["outdir"]+"/{prefix}/QC_workflow/gunc/{sample}.temp"))
	conda:
		"envs/gunc.yaml"
	log:
		config["base_log_outdir"]+"/{prefix}/QC_workflow/gunc/logs/{sample}/gunc.log"
	shell:
		"""
		mkdir {output.tempdir}
		gunc run -i {input.assembly} --file_suffix .fasta --use_species_level --temp_dir {output.tempdir} --detailed_output --out_dir {output.guncoutdir} -r gunc_db/gunc_db_progenomes2.1.dmnd
		"""


rule run_gunc_summarise:
	input:
		#gunc report files
		gunc_reports=expand(config["outdir"]+"/{prefix}/QC_workflow/gunc/{sample}.out/GUNC.progenomes_2.1.maxCSS_level.tsv", prefix=prefix, sample=sample_ids)
	output:
		#combined, tab-delimited gunc data
		combined_gunc=config["outdir"]+"/{prefix}/QC_workflow/summaries/gunc_report.txt"
	log:
		"logs/{prefix}/summaries/combine_gunc.log"
	threads: 1
	script:
		"../../scripts/combine_gunc.py"