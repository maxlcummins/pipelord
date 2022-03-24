import re
import subprocess
import os

prefix = config["prefix"]
outdir = config["outdir"]
maxthreads = snakemake.utils.available_cpu_count()

if path.exists(config["gunc_db_path"]) == False:
	if path.exists('resources/dbs/gunc_db') == False:
		print("Gunc DB not detected. It will be downloaded (to 'resources/dbs/gunc_db'), which will take some time as it is (~13GB)")
		rule gunc_db_dl:
			output:
				gunc_db = directory("resources/dbs/gunc_db"),
			conda:
				"../envs/gunc.yaml"
			shell:
				"""
				mkdir -p resources/dbs/gunc_db
				gunc download_db {output}
				"""
	gunc_db_path = 'resources/dbs/gunc_db'
else: gunc_db_path = config["gunc_db_path"]

rule gunc_run:
	input:
		gunc_db = gunc_db_path,
		assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
	output:
		guncoutdir = directory(config["outdir"]+"/{prefix}/QC_workflow/gunc/{sample}.out"),
		tempdir = directory(temp(config["outdir"]+"/{prefix}/QC_workflow/gunc/{sample}.temp"))
	conda:
		"../envs/gunc.yaml"
	log:
		config["base_log_outdir"]+"/{prefix}/QC_workflow/gunc/logs/{sample}/gunc.log"
	shell:
		"""
		mkdir -p {output.guncoutdir}
		mkdir -p {output.tempdir}
		gunc run -i {input.assembly} --file_suffix .fasta --use_species_level --temp_dir {output.tempdir} --detailed_output --out_dir {output.guncoutdir} -r {input.gunc_db}/gunc_db_progenomes2.1.dmnd
		"""


rule run_gunc_summarise:
	input:
		#gunc report files
		gunc_reports=expand(config["outdir"]+"/{prefix}/QC_workflow/gunc/{sample}.out", prefix=prefix, sample=sample_ids)
	output:
		#combined, tab-delimited gunc data
		combined_gunc=config["outdir"]+"/{prefix}/QC_workflow/summaries/gunc_report.txt"
	log:
		"logs/{prefix}/summaries/combine_gunc.log"
	threads: 1
	script:
		"../../scripts/combine_gunc.py"
