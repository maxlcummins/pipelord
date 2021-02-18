import re
import subprocess
import os

#configfile: "misc/masterconfig2.yaml"

#print(assembly_path)

# Get reads
#sample_ids, = glob_wildcards(config['raw_reads_path']+"/{sample}.R1.fastq.gz")
prefix = config['prefix']
maxthreads = snakemake.utils.available_cpu_count()

#print(config['krakendb'])
#print(sample_ids)

# one rule to rule them all
#rule all:
#    input:
#        expand(config['outdir']+"/{prefix}/kraken2/{sample}.out", sample=sample_ids, prefix=prefix),
#        expand(config['outdir']+"/{prefix}/kraken2/{sample}.report",sample=sample_ids, prefix=prefix),
#        expand(config['outdir']+"/{prefix}/summaries/kraken2_summary.txt", prefix=prefix),


rule run_kraken2:
    input:
        db = config['krakendb'],
        assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        out = temp(config['outdir']+"/{prefix}/kraken2/{sample}.out"),
        report = temp(config['outdir']+"/{prefix}/kraken2/{sample}.report")
    log:
        config['base_log_outdir']+"/{prefix}/kraken2/run/{sample}.log"
    conda:
        "config/kraken2.yaml"
    shell:
        """
        kraken2 --db {input.db} --use-names --report {output.report} --output {output.out} {input.assembly}  2>&1 {log}
        """

rule name_kraken2_inline:
    input:
        config['outdir']+"/{prefix}/kraken2/{sample}.report"
    output:
        temp(config['outdir']+"/{prefix}/kraken2/{sample}_clean.report")
    shell:
        """
        awk 'NR == 1 {{print FILENAME "\t" $0; next;}}{{print FILENAME "\t" $0;}}' {input} > {output}
        perl -p -i -e 's@\.report@@g' {output}
        """

rule combine_kraken2_reports:
    input:
        expand(config['outdir']+"/{prefix}/kraken2/{sample}_clean.report",sample=sample_ids, prefix=prefix)
    output:
        combined_output = temp(config['outdir']+"/{prefix}/summaries/kraken2_summary_temp.txt"),
        combined_output_cleaned = temp(config['outdir']+"/{prefix}/summaries/kraken2_summary_temp2.txt"),
        simplified_combined = temp(config['outdir']+"/{prefix}/summaries/kraken2_summary_temp3.txt"),
        summary_full = config['outdir']+"/{prefix}/summaries/kraken2_full_summary.txt",
        summary = config['outdir']+"/{prefix}/summaries/kraken2_summary.txt"

    shell:
        """
        # Combine the summaries
        cat {input} > {output.combined_output}
        # Clean sample names
        perl -p -i -e s'@.*kraken2/@@g' {output.combined_output}
        # Removes duplicate headers (in this case lines starting with filename)
        awk 'FNR==1 {{ header = $0; print }} $0 != header' {output.combined_output} > {output.combined_output_cleaned}
        # Create simple summary - 1. find cols with species designations (== S) 2. Sort by highest match in column 2. 3. Select the row with the first instance of each sample (i.e. highest species match)
        grep -P "\tS\t" {output.combined_output_cleaned} | sort -nk2 -r | awk '{{ if (a[$1]++ == 0) print $0; }}' $@ > {output.simplified_combined}
        #awk '{{if ($5 == "S") print $0;}}' {output.combined_output_cleaned} | sort -nk2 -r | awk '{{ if (a[$1]++ == 0) print $0; }}' $@ > {output.simplified_combined}
        # Insert a header
        echo -e "name\tperc_frags\tnum_frags\tnum_frags_direct\trank\tncbi_taxon_id\tscientific_name" | cat - {output.simplified_combined} > {output.summary}
        # Insert a header
        echo -e "name\tperc_frags\tnum_frags\tnum_frags_direct\trank\tncbi_taxon_id\tscientific_name" | cat - {output.combined_output_cleaned} > {output.summary_full}
        """

#include: "genome_assembly.smk"