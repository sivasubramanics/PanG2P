
def get_samples_map():
    samples = {}
    with open(config["samples"], "r") as f:
        # skip header
        next(f)
        for line in f:
            sample_id, path = line.strip().split(",")
            samples[sample_id] = path
    return samples


SAMPLES = get_samples_map()

def get_fq_dir(wildcards):
    return SAMPLES[wildcards.sample]

rule bbduk_kmer_search:
    input:
        ref="work/03.pankmers/Pan_kmers.fasta",
        query=get_fq_dir
    output:
        temp("work/04.genotyping/{sample}/bbduk.tsv")
    conda:
        ENV_DIR / "genotyping.yaml"
    container:
        CONTAINER
    log:
        "logs/04.genotyping/bbduk_kmer_search/{sample}.log"
    params:
        kmer_size = config["pan_kmer_size"],
        memory = config["max_mem"]
    benchmark:
        "benchmarks/04.genotyping/bbduk_kmer_search/{sample}.tsv"
    shell:
        r"""
            (
                mkdir -p $(dirname {output})

                find {input.query} | egrep 'fq.gz$|fastq.gz$|R1.fastq.gz$' \
                | xargs gzip -cd \
                | bbduk.sh \
                    -Xmx{params.memory}G \
                    -k={params.kmer_size} \
                    overwrite=t \
                    int=f \
                    ref={input.ref} \
                    in=stdin.fq \
                    out=/dev/null \
                    rpkm={output} \
                    nzo=f
            ) > {log} 2>&1
        """

rule summarize_kmer_counts:
    input:
        "work/04.genotyping/{sample}/bbduk.tsv"
    output:
        "work/04.genotyping/{sample}/count.tsv.gz"
    log:
        "logs/04.genotyping/summarize_kmer_counts/{sample}.log"
    benchmark:
        "benchmarks/04.genotyping/summarize_kmer_counts/{sample}.tsv"
    shell:
        r"""
        awk '
        BEGIN {{ OFS="\t" }}

        /^#/ {{ next }}

        {{
            if ($5 < 1)
                next

            key = $1
            sub(/::.*/, "", key)
            sub(/_[0-9]+$/, "", key)

            n[key]++
            sum[key] += $5
        }}

        END {{
            for (k in n)
                print k, n[k], sum[k]
        }}
        ' {input} \
        | sort \
        | gzip > {output}
        """

rule collect_kmer_counts:
    input:
        expand(
            "work/04.genotyping/{sample}/count.tsv.gz",
            sample=SAMPLES.keys()
        ),
        "work/03.pankmers/Pan_count_kmers.tsv"
    output:
        "work/04.genotyping/overall_counts.tsv.gz",
        "work/04.genotyping/overall_mapped_kmers.tsv.gz"
    params:
        script = SCRIPTS_DIR / "collect_kmer_counts.R",
        input_dir = "work/04.genotyping/",
        output_dir = "work/04.genotyping/",
    threads:
        config["threads"]
    conda:
        ENV_DIR / "genotyping.yaml"
    container:
        CONTAINER
    log:
        "logs/04.genotyping/collect_kmer_counts.log"
    benchmark:
        "benchmarks/04.genotyping/collect_kmer_counts.tsv"
    shell:
        r"""
        Rscript {params.script} -i {params.input_dir} -p {input[1]} -o {params.output_dir} -c {threads} > {log} 2>&1
        """

rule normalize_kmer_counts:
    input:
        overall = "work/04.genotyping/overall_counts.tsv.gz",
        mapped = "work/04.genotyping/overall_mapped_kmers.tsv.gz"
    output:
        "work/04.genotyping/genotype_matrix.tsv.gz",
        "work/04.genotyping/genotype_map.tsv.gz"
    params:
        script = SCRIPTS_DIR / "normalize_kmer_counts.R",
        input_dir = "work/04.genotyping/",
        output_dir = "work/04.genotyping/"
    conda:
        ENV_DIR / "genotyping.yaml"
    container:
        CONTAINER
    log:
        "logs/04.genotyping/normalize_kmer_counts.log"
    benchmark:
        "benchmarks/04.genotyping/normalize_kmer_counts.tsv"
    shell:
        r"""
        Rscript {params.script} --cm {input.mapped} --oc {input.overall} --outdir {params.output_dir} > {log} 2>&1
        """

rule genotyping_complete:
    input:
        "work/04.genotyping/genotype_matrix.tsv.gz",
        "work/04.genotyping/genotype_map.tsv.gz"
    output:
        "work/04.genotyping/.done"
    log:
        "logs/04.genotyping/genotyping_complete.log"
    benchmark:
        "benchmarks/04.genotyping/genotyping_complete.tsv"
    shell:
        "touch {output}"
