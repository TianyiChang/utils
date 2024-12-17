mge_aa_split_batch, = glob_wildcards("../query_seqs/{mge_aa_split_batch}.fasta")

rule all:
    input:
        expand(
            "../pfama_hmm_clean/{mge_aa_split_batch}.done",
            mge_aa_split_batch=mge_aa_split_batch
            )

rule hmmsearch_pfama:
    input:
        query="../query_seqs/{mge_aa_split_batch}.fasta",
        ref=config["path2refdb"] + "/pfam_a/Pfam-A.hmm"
    output:
        touch("../pfama_hmm/{mge_aa_split_batch}.done")
    params:
        bitscore=10,
        output="../pfama_hmm/{mge_aa_split_batch}.tblout"
    threads: 2
    resources:
        mem_mb=5120
    shell:
        """
        source {config[path2conda]}
        conda activate {config[path2scripts]}/../conda/hmmer_pfam_annot_hmmer

        hmmsearch --cpu {threads} -T {params.bitscore} \
            --tblout {params.output} \
            {input.ref} {input.query}
        """

rule clean_hmm_out:
    input:
        checkpoint="../pfama_hmm/{mge_aa_split_batch}.done"
    output:
        touch("../pfama_hmm_clean/{mge_aa_split_batch}.done")
    threads: 1
    params:
        pfam_in="../pfama_hmm/{mge_aa_split_batch}.tblout",
        pfam_out="../pfama_hmm_clean/{mge_aa_split_batch}.txt"
    resources:
        mem_mb=1024
    shell:
        """
        echo "gene pfam_name pfam_id evalue score" \
            > {params.pfam_out}

        grep -v "^#" {params.pfam_in} | \
            awk '{{print $1, $3, $4, $5, $6}}' >> {params.pfam_out}

        rm {params.pfam_in}
        """
        