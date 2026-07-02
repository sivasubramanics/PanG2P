rule mmseqs_createdb:
    input:
        "work/01.single_copy/single_copy.fa"
    output:
        "work/02.clustering/single_copy.db",
        temp(directory("work/02.clustering/tmp"))
    conda:
        ENV_DIR / "clustering.yaml"
    container:
        CONTAINER
    threads: config["threads"]
    log:
        "logs/02.clustering/mmseqs_createdb.log"
    benchmark:
        "benchmarks/02.clustering/mmseqs_createdb.tsv"
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
        "work/02.clustering/tmp"
    output:
        "work/02.clustering/clusters.tsv"
    conda:
        ENV_DIR / "clustering.yaml"
    container:
        CONTAINER
    params:
        prefix = "work/02.clustering/clusters",
        min_seq_id = config["min_seq_id"],
        min_coverage = config["min_coverage"],
        coverage_mode = config["coverage_mode"]
    threads: config["threads"]
    log:
        "logs/02.clustering/mmseqs_clustering.log"
    benchmark:
        "benchmarks/02.clustering/mmseqs_clustering.tsv"
    shell:
        r"""
            (
                rm -rf {input[1]} {params.prefix}.*
                mmseqs cluster --threads {threads} {input[0]} {params.prefix} {input[1]} --min-seq-id {params.min_seq_id} -c {params.min_coverage} --cov-mode {params.coverage_mode}
                mmseqs createtsv --threads {threads} {input[0]} {input[0]} {params.prefix} {output}
            ) > {log} 2>&1
        """

rule get_cluster_representatives:
    input:
        clusters="work/02.clustering/clusters.tsv",
        fai="work/01.single_copy/single_copy.fa.fai"
    output:
        members="work/02.clustering/pan_cluster_all_member.txt",
        reps="work/02.clustering/pan_cluster_mmseq_rep.txt",
        bed="work/02.clustering/pan_cluster_bed_rep.txt"
    conda:
        ENV_DIR / "r.yaml"
    container:
        CONTAINER
    params:
        script=SCRIPT_DIR / "cluster_representatives.R"
    log:
        "logs/02.clustering/get_cluster_representatives.log"
    benchmark:
        "benchmarks/02.clustering/get_cluster_representatives.tsv"
    shell:
        r"""
        Rscript {params.script} \
            --clusters {input.clusters} \
            --fai {input.fai} \
            --outdir $(dirname {output.members})
        """

rule clustering_complete:
    input:
        members="work/02.clustering/pan_cluster_all_member.txt",
        reps="work/02.clustering/pan_cluster_mmseq_rep.txt",
        bed="work/02.clustering/pan_cluster_bed_rep.txt"
    output:
        "work/02.clustering/.done"
    log:
        "logs/02.clustering/clustering_complete.log"
    benchmark:
        "benchmarks/02.clustering/clustering_complete.tsv"
    shell:
        "touch {output}"
