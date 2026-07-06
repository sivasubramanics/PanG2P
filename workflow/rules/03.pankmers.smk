def get_genomes_list(wildcards):
    return list(genome_map.keys())

rule split_and_sort_pan_clusters:
    input:
        "work/02.clustering/pan_cluster_bed_rep.txt"
    output:
        "work/03.pankmers/{genome}.tsv"
    log:
        "logs/03.pankmers/split_and_sort_pan_clusters/{genome}.log"
    benchmark:
        "benchmarks/03.pankmers/split_and_sort_pan_clusters/{genome}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output})

        tail -n +2 {input} |
        grep -v 'scaffold' |
        awk -F'\t' -v grp="{wildcards.genome}" '$5 == grp' |
        sort -k1,1 -k2,2n \
        > {output}
        """

rule make_intervals:
    input:
        "work/03.pankmers/{genome}.tsv"
    output:
        "work/03.pankmers/{genome}_interval_formatted.bed"
    conda:
        ENV_DIR / "pankmers.yaml"
    container:
        CONTAINER
    params:
        script=SCRIPTS_DIR / "reformat_kmer.py",
        kmer_size=config["pan_kmer_size"],
        step_size=config["pan_kmer_step"]
    log:
        "logs/03.pankmers/make_intervals/{genome}.log"
    benchmark:
        "benchmarks/03.pankmers/make_intervals/{genome}.tsv"
    shell:
        r"""
        tmp=$(mktemp)

        bedtools makewindows -b {input} -w {params.kmer_size} -s {params.step_size} -i src \
        | awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4,($3-$2)}}' \
        | awk '$5=={params.kmer_size}' \
        | cut -f1-4 > "$tmp"

        python {params.script} "$tmp" {output}

        rm -f "$tmp"
        """

rule extract_kmers:
    input:
        fasta="work/01.single_copy/{genome}.fasta",
        bed="work/03.pankmers/{genome}_interval_formatted.bed"
    output:
        "work/03.pankmers/{genome}_sc_kmer.fasta"
    conda:
        ENV_DIR / "pankmers.yaml"
    container:
        CONTAINER
    log:
        "logs/03.pankmers/extract_kmers/{genome}.log"
    benchmark:
        "benchmarks/03.pankmers/extract_kmers/{genome}.tsv"
    shell:
        r"""
        cut -f1-4 {input.bed} \
        | bedtools getfasta \
            -name \
            -fi {input.fasta} \
            -bed /dev/stdin \
            -fo {output}
        """

rule combine_kmers:
    input:
        expand(
            "work/03.pankmers/{genome}_sc_kmer.fasta",
            genome=get_genomes_list
        )
    output:
        "work/03.pankmers/Pan_kmers.fasta"
    log:
        "logs/03.pankmers/combine_kmers.log"
    benchmark:
        "benchmarks/03.pankmers/combine_kmers.tsv"
    shell:
        """
        cat {input} > {output}
        """

rule count_kmers:
    input:
        "work/03.pankmers/Pan_kmers.fasta"
    output:
        "work/03.pankmers/Pan_count_kmers.tsv"
    log:
        "logs/03.pankmers/count_kmers.log"
    benchmark:
        "benchmarks/03.pankmers/count_kmers.tsv"
    shell:
        r"""
        awk '
        /^>/ {{
            key = substr($1, 2)
            sub(/::.*/, "", key)
            sub(/_[0-9]+$/, "", key)
            count[key]++
        }}

        END {{
            for (k in count)
                print k "\t" count[k]
        }}
        ' {input} \
        | sort \
        > {output}
        """


rule complete_pankmers:
    input:
        "work/03.pankmers/Pan_kmers.fasta",
        "work/03.pankmers/Pan_count_kmers.tsv"
    output:
        "work/03.pankmers/.done"
    log:
        "logs/03.pankmers/complete_pankmers.log"
    benchmark:
        "benchmarks/03.pankmers/complete_pankmers.tsv"
    shell:
        r"""
        touch {output}
        """
