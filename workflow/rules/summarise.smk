#rule summarise_all:
#    input:
#        mlst = expand(config['outdir']+"/{prefix}/summaries/mlst.txt", prefix=prefix),
#        amr_snps = expand(config['outdir']+"/{prefix}/summaries/Pointfinder.txt", prefix=prefix),
#        genotype = expand(config['outdir']+"/{prefix}/summaries/abricate_hits.txt", prefix=prefix),
#        pmlst = expand(config['outdir']+"/{prefix}/summaries/pMLST.txt", prefix=prefix),
#        #species_id = expand(config['outdir']+"/{prefix}/summaries/kraken2_full_summary.txt", prefix = prefix),
#        #assembly_stats = expand(config['outdir']+"/{prefix}/summaries/assembly_stats.txt", prefix = prefix),
#    output:
#        expand(config['outdir']+"/{prefix}/summaries/{prefix}_simple_summary_N"+str(config['abricateR_identity'])+"L"+str(config['abricateR_length'])+".csv", prefix=prefix)
#    conda:
#        "../envs/R.yaml"
#    log:
##        out = config['base_log_outdir']+"/{prefix}/summarise/summarise_out.log",
##        err = config['base_log_outdir']+"/{prefix}/summarise/summarise_err.log"
#    params:
#        prefix = config['prefix'],
#        identity = config['abricateR_identity'],
#        length = config['abricateR_length'],
#        outdir = expand(config['outdir']+"/{prefix}/summaries/", prefix=prefix),
#        email = config['email']
#    shell:
#        """
#        Rscript scripts/abricateR.R '{input.genotype}' '{params.prefix}' '{params.outdir}' '{params.identity}' '{params.length}' '{input.amr_snps}' '{input.pmlst}'
#        echo Your Snakemake job with prefix \'{params.prefix}\' has finished running. | mail -s 'Snakemake job has finished' {params.email}
#        """