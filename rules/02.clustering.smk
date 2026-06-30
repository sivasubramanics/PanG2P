rule mmseqs_createdb:
    input:
        "work/01.single_copy/single_copy.fa"
    output:
        "work/02.clustering/single_copy.db",
        temp(directory("work/02.clustering/tmp"))
    conda:
        "../envs/clustering.yaml"
    threads: config["threads"]
    log:
        "logs/02.clustering/single_copy.createdb.log"
    shell:
        r"""
            (
                mkdir -p {output[1]}
                mmseqs createdb {input} {output[0]} --threads {threads} > {log} 2>&1
            ) > {log} 2>&1
        """

rule mmseqs_clustering:
    input:
        "work/02.clustering/single_copy.db",
        directory("work/02.clustering/tmp")
    output:
        "work/02.clustering/clusters.tsv"
    conda:
        "../envs/clustering.yaml"
    params:
        prefix = "work/02.clustering/clusters",
        min_seq_id = config["min_seq_id"],
        min_coverage = config["min_coverage"],
        coverage_mode = config["coverage_mode"]
    threads: config["threads"]
    log:
        "logs/02.clustering/clusters.log"
    shell:
        r"""
            (
                mmseqs cluster --threads {threads} {input[0]} {params.prefix} {input[1]} --min-seq-id {params.min_seq_id} -c {params.min_coverage} --cov-mode {params.coverage_mode}
                mmseqs createtsv --threads {threads} {params.prefix} {params.prefix} {params.prefix} {output}
            ) > {log} 2>&1
        """

rule clustering_complete:
    input:
        "work/02.clustering/clusters.tsv"
    output:
        "work/02.clustering/.done"
    shell:
        "touch {output}"
