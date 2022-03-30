rule roary_subset:
    input:
        config['outdir']+"/{prefix}/dfast/gffs/{subset}.gff"
    output:
        config['outdir']+"/{prefix}/dfast/{subset_name}/{subset}.gff"
    shell:
        """
        cp -n {input} {output}
        """

rule roary:
    input:
        expand(config['outdir']+"/{prefix}/dfast/{subset_name}/{subset}.gff", prefix=prefix, subset_name=subset_name, subset=subset)
    output:
        directory(config['outdir']+"/{prefix}/roary/{subset_name}")
    conda:
        "../envs/roary.yaml"
    log:
        config['base_log_outdir']+"/{prefix}/roary/{subset_name}.log"
    threads: 12
    params:
        out = config['outdir']+"/{prefix}/roary/{subset_name}",
        core_req = config['roary_core_req'],
        kraken_db = config['krakendb']
    shell:
        "roary -p {threads} -e -v -n -r -cd {params.core_req} -qc -k {params.kraken_db} -f {params.out} {input}  2> {log}"

rule snp_sites:
    input:
        config['outdir']+"/{prefix}/roary/{subset_name}/core_gene_alignment.aln"
    output:
        config['outdir']+"/{prefix}/roary/{subset_name}/core_gene_alignment_snp_sites.aln"
    conda:
        "config/snp_sites.yaml"
    threads:
        1
    log:
        config['base_log_outdir']+"/{prefix}/roary/logs/{subset_name}/snp_sites.log"
    shell:
        "snp-sites -c {input} > {output}"

rule iqtree:
    input:
        config['outdir']+"/{prefix}/roary/{subset_name}/core_gene_alignment.aln"
    output:
        config['outdir']+"/{prefix}/trees/{subset_name}/roary_"+str(config['roary_core_req'])+"/core_gene_alignment.aln.treefile"
    conda:
        "config/iqtree.yaml"
    threads:
        12
    log:
        config['base_log_outdir']+"/{prefix}/roary/logs/{subset_name}/iqtree.log"
    shell:
        "iqtree -s {input} -m MFP -ntmax -bb 1000"


rule iqtree_snp_sites:
    input:
        config['outdir']+"/{prefix}/roary/{subset_name}/core_gene_alignment_snp_sites.aln"
    output:
        config['outdir']+"/{prefix}/pangenome/{subset_name}/roary_"+str(config['roary_core_req'])+"/core_gene_alignment_snp_sites.aln.treefile"
    conda:
        "config/iqtree.yaml"
    threads:
        12
    log:
        config['base_log_outdir']+"/{prefix}/roary/logs/{subset_name}/iqtree_snp_sites.log"
    shell:
        "iqtree -s {input} -m MFP -ntmax -bb 1000"

#snp-dists -c full_aln/core_gene_alignment.aln > full_aln/core_gene_alignment.csv
#snp-dists -c snp_sites/core_gene_alignment_snp_sites.aln > snp_sites/core_gene_alignment_snp_sites.csv