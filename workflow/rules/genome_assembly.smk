prefix = config['prefix']

logs = config['base_log_outdir']

if config['input_type'] ==  "raw_reads":
    assembly_path = config['outdir']+"/"+config['prefix']+"/shovill/assemblies"
else:
    assembly_path = config['assembly_path']


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
            config['outdir']+"/{prefix}/shovill/assembly_stats/{sample}_assembly_stats.txt"
        conda:
            "../envs/assembly_stats.yaml"
        shell:
            "assembly-stats -t {input} > {output}"


elif config['input_type'] == 'raw_reads':
    rule run_shovill:
        input:
            r1_filt = config['outdir']+"/{prefix}/fastp/{sample}.R1.fastq.gz",
            r2_filt = config['outdir']+"/{prefix}/fastp/{sample}.R2.fastq.gz"
        output:
            shov_out = directory(config['outdir']+"/{prefix}/shovill/shovill_out/{sample}.out"),
            assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
        log:
            out = config['base_log_outdir']+"/{prefix}/shovill/run/{sample}_out.log",
            err = config['base_log_outdir']+"/{prefix}/shovill/run/{sample}_err.log"
        threads:
            8
        conda:
            "../envs/shovill.yaml"
        shell:
            """
            shovill --minlen 200 --outdir {output.shov_out} --R1 {input.r1_filt} --R2 {input.r2_filt} 1> {log.out} 2> {log.err}
            cp {output.shov_out}/contigs.fa {output.assembly}
            """
    rule run_assembly_stats:
        input:
            config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
        output:
            config['outdir']+"/{prefix}/shovill/assembly_stats/{sample}_assembly_stats.txt"
        conda:
            "../envs/assembly_stats.yaml"
        shell:
            "assembly-stats -t {input} > {output}"
