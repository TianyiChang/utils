rule all:
    input:
        config["path2outputs"] + "/checkpoints/split_combined_mge_aa_into_chunks.done"

rule split_combined_4_annot:
    input:
        pyscript=config["path2scripts"] + "/split_fa_by_chars.py",
        smkscript1=config["path2scripts"] + "/mge_protein_annot.smk"
    output:
        checkpoint=touch(config["path2outputs"] +
            "/checkpoints/split_combined_4_annot.done")
    threads: 1
    params:
        infile=config["path2outputs"] + "/checkpoints/query.faa",
        outdir=config["path2outputs"] + "/checkpoints/candid_mge_aa_chunks",
        chars_per_file=3_000_000,
        files_per_folder=1_000
    resources:
        mem_mb=20480
    shell:
        """
        source {config[path2conda]}
        conda activate {config[path2scripts]}/../conda/hmmer_pfam_annot_python3
        rm -r {params.outdir} || true
        mkdir -p {params.outdir}        

        python {input.pyscript} \
            -i {params.infile} -o {params.outdir} \
            --scripts {input.smkscript1} \
            --filesize {params.chars_per_file} \
            --foldersize {params.files_per_folder}
        """

rule remove_query:
    input:
        config["path2outputs"] + "/checkpoints/split_combined_4_annot.done"
    output:
        checkpoint=touch(config["path2outputs"] +
            "/checkpoints/split_combined_mge_aa_into_chunks.done")
    threads: 1
    params:
        seq=config["path2outputs"] + "/checkpoints/query.faa",
    resources:
        mem_mb=10
    shell:
        """
        rm {params.seq}
        """
