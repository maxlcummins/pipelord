prefix = config['prefix']

logs = config['base_log_outdir']

if config['input_type'] == "assemblies":
    rule pseudo_assemble:
        input:
            assembly = config['genome_path']+"/{sample}.fasta"
        output:
            assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
        threads:
            1
        shell:
            """
            cp -n {input} {output}
            """

elif config['input_type'] == "reads":
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
            shovill --minlen 200 --cpus {threads} --outdir {output.shov_out} --R1 {input.r1_filt} --R2 {input.r2_filt} 1> {log.out} 2> {log.err}
            cp {output.shov_out}/contigs.fa {output.assembly}
            """
else:
    print("Error: input_type must be set to either 'reads' or 'assemblies'. Please edit the config file accordingly.")

rule run_assembly_stats:
    input:
        config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
    output:
        config['outdir']+"/{prefix}/shovill/assembly_stats/{sample}_assembly_stats.txt"
    conda:
        "../envs/assembly_stats.yaml"
    threads:
        2
    shell:
        "assembly-stats -t {input} > {output}"

rule touch_assemblies: #required for rules where we need a directory of fastas as an input
    input:
        expand(config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta", prefix=prefix, sample=sample_ids)
    output:
        assemblydir = directory(config['outdir']+"/{prefix}/shovill/assemblies_temp")
    threads: 1
    shell:
        """
        mkdir {output}
        cp {input} {output}
        """

rule run_assembly_summarise:
    input:
        assembly_stats_summary=expand(config['outdir']+"/{prefix}/shovill/assembly_stats/{sample}_assembly_stats.txt", prefix=prefix, sample=sample_ids)
    output:
        combine_assembly_stats=config["outdir"]+"/{prefix}/summaries/assembly_stats.txt"
    params:
        extra="",
    log:
        "logs/{prefix}/summaries/combine_assembly_stats.log",
    threads: 1
    script:
        "../../scripts/combine_assembly_stats.py"
