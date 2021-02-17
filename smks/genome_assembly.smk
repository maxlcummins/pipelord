#configfile: "misc/masterconfig2.yaml"

#sample_ids, = glob_wildcards(config['raw_reads_path']+"/{sample}.R1.fastq.gz")

#print(sample_ids)
prefix = config['prefix']

logs = config['base_log_outdir']

if config['input_type'] ==  "raw_reads":
    assembly_path = config['outdir']+"/"+config['prefix']+"/shovill/assemblies"
else:
    assembly_path = config['assembly_path']

#rule all:
#    input:
#        expand(config['outdir']+"/{prefix}/shovill/shovill_out/{sample}.out", sample=sample_ids, prefix=prefix),
#        expand(config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta", sample=sample_ids, prefix=prefix),
#        expand(config['outdir']+"/{prefix}/shovill/assembly_stats/{sample}_assembly_stats.txt", sample=sample_ids, prefix=prefix),
#        #expand(config['outdir']+"/{prefix}/summaries/assembly_stats.txt", prefix=prefix)


if config['input_type'] == 'assemblies':
    rule pseudo_assemble:
        input:
            assembly = config['assembly_path']+"/{sample}.fasta"
        output:
            assembly = temp(config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta")
        shell:
            "cp {input} {output}"

    rule run_assembly_stats_:
        input:
            config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
        output:
            temp(config['outdir']+"/{prefix}/shovill/assembly_stats/{sample}_assembly_stats.txt")
        conda:
            "config/assembly_stats.yaml"
        shell:
            "assembly-stats -t {input} > {output}"

    rule combine_assembly_stats_:
        input:
            expand(config['outdir']+"/{prefix}/shovill/assembly_stats/{sample}_assembly_stats.txt", sample=sample_ids, prefix=prefix)
        output:
            stats_temp = temp(config['outdir']+"/{prefix}/summaries/temp_assembly_stats.txt")
        shell:
            """
            cat {input} > {output}
            """

    rule clean_assembly_stats_:
        input:
            stats_temp = config['outdir']+"/{prefix}/summaries/temp_assembly_stats.txt",
        output:
            stats = config['outdir']+"/{prefix}/summaries/assembly_stats.txt"
        shell:
            """
            # Cleans file names
            perl -p -i -e 's@.*assemblies/@@g' {input}
            perl -p -i -e 's@.fasta@@g' {input}
            # Changes column name to name rather than filename
            perl -p -i -e 's@^filename@name@g' {input}
            #Removes duplicate headers (in this case lines starting with filename)
            awk 'FNR==1 {{ header = $0; print }} $0 != header' {input} > {output}
            """

elif config['input_type'] == 'raw_reads':
    rule run_shovill:
        input:
            r1_filt = config['outdir']+"/{prefix}/fastp/{sample}.R1.fastq.gz",
            r2_filt = config['outdir']+"/{prefix}/fastp/{sample}.R2.fastq.gz"
        output:
            shov_out = temp(directory(config['outdir']+"/{prefix}/shovill/shovill_out/{sample}.out")),
            assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
        log:
            config['base_log_outdir']+"/{prefix}/shovill/run/{sample}_out.log"
        threads:
            4
        conda:
            "config/shovill.yaml"
        shell:
            """
            shovill --minlen 200 --outdir {output.shov_out} --R1 {input.r1_filt} --R2 {input.r2_filt} 1> {log}
            cp {output.shov_out}/contigs.fa {output.assembly}
            """
    rule run_assembly_stats:
        input:
            config['outdir']+"/{prefix}/shovill/shovill_out/{sample}.out"
        output:
            temp(config['outdir']+"/{prefix}/shovill/assembly_stats/{sample}_assembly_stats.txt")
        conda:
            "config/assembly_stats.yaml"
        shell:
            "assembly-stats -t {input}/contigs.fa > {output}"

    rule combine_assembly_stats:
        input:
            expand(config['outdir']+"/{prefix}/shovill/assembly_stats/{sample}_assembly_stats.txt", sample=sample_ids, prefix=prefix)
        output:
            stats_temp = temp(config['outdir']+"/{prefix}/summaries/temp_assembly_stats.txt")
        shell:
            """
            cat {input} > {output}
            """

    rule clean_assembly_stats:
        input:
            stats_temp = config['outdir']+"/{prefix}/summaries/temp_assembly_stats.txt",
        output:
            stats = config['outdir']+"/{prefix}/summaries/assembly_stats.txt"
        shell:
            """
            # Cleans file names
            perl -p -i -e 's@.*shovill_out/@@g' {input}
            perl -p -i -e 's@.out/contigs.fa@@g' {input}
            # Changes column name to name rather than filename
            perl -p -i -e 's@^filename@name@g' {input}
            #Removes duplicate headers (in this case lines starting with filename)
            awk 'FNR==1 {{ header = $0; print }} $0 != header' {input} > {output}
            """

#include: "read_cleaning.smk"