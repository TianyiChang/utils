rule all:
    input:
        config["path2outputs"] +
            "/checkpoints/join_annot.done"

rule join_annot:
    input:
        rscript=config["path2scripts"] + "/join_annot.r"
    output:
        touch(config["path2outputs"] +
            "/checkpoints/join_annot.done")
    threads: 30
    params:
        checkpdir=config["path2outputs"] + "/checkpoints",
        pfam_metad=config["path2refdb"] + "/pfam_a/Pfam-A.clans.tsv",
        outfile=config["path2outputs"] + "/pfam_annot.tsv"
    resources:
        mem_mb=20480
    shell:
        """
        source {config[path2conda]}
        conda activate {config[path2scripts]}/../conda/hmmer_pfam_annot_R_pcks

        Rscript {input.rscript} {params.checkpdir} {params.pfam_metad} {params.outfile}
        """
