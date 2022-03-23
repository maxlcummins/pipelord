prefix = config['prefix']

logs = config['base_log_outdir']

if config["genotype_modules"]["run_genome_assembly"] == False:
    rule pseudo_assemble:
        input:
            assembly = config['genome_path']+"/{sample}.fasta"
        output:
            assembly = config['outdir']+"/{prefix}/shovill/assemblies/{sample}.fasta"
        shell:
            """
            cp -n {input} {output}
            """

elif config["genotype_modules"]["run_fastp"] == False:
    rule run_shovill:
        input:
            r1 = config['genome_path']+"/{sample}.R1.fastq.gz",
            r2 = config['genome_path']+"/{sample}.R2.fastq.gz"
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
            shovill --minlen 200 --outdir {output.shov_out} --R1 {input.r1} --R2 {input.r2} 1> {log.out} 2> {log.err}
            cp {output.shov_out}/contigs.fa {output.assembly}
            """
else:
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
