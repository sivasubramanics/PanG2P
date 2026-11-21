
def get_samples_map():
    samples = {}
    with open(config["samples"], "r") as f:
        next(f)  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            sample_id, path = line.split(",", 1)
            samples.setdefault(sample_id, []).append(path)
    return samples


SAMPLES = get_samples_map()


def get_fq_files(wildcards):
    return SAMPLES[wildcards.sample]

rule bbduk_kmer_search:
    input:
        ref="work/03.pankmers/Pan_kmers.fasta",
        query=get_fq_files
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
    resources:
        memory = int(config["max_mem"] * 0.9)
    threads:
        config["threads"]
    benchmark:
        "benchmarks/04.genotyping/bbduk_kmer_search/{sample}.tsv"
    shell:
        r"""(
            set -euo pipefail
            mkdir -p "$(dirname "{output}")"

            files=()
            for f in {input.query}; do
                case "$f" in
                    *.fq.gz|*.fastq.gz)
                        if [ ! -f "$f" ]; then
                            echo "ERROR: Missing input FASTQ inside container: $f" >&2
                            exit 1
                        fi
                        files+=("$f")
                        ;;
                    *)
                        echo "ERROR: Unsupported file extension: $f" >&2
                        exit 1
                        ;;
                esac
            done

            if [ "${{#files[@]}}" -eq 0 ]; then
                echo "ERROR: No FASTQ files provided for {wildcards.sample}" >&2
                exit 1
            fi

            gzip -ckd "${{files[@]}}" \
            | bbduk.sh \
                -Xmx{resources.memory}G \
                -k={params.kmer_size} \
                overwrite=t \
                int=f \
                ref={input.ref} \
                in=stdin.fq \
                out=/dev/null \
                rpkm={output} \
                nzo=f \
                t={threads}
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
            sub(/_L$/, "", key)

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
        count_files = expand(
            "work/04.genotyping/{sample}/count.tsv.gz",
            sample=SAMPLES.keys()
        ),
        pan_kmer_counts = "work/03.pankmers/Pan_count_kmers.tsv"
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
        Rscript {params.script} -i {params.input_dir} -p {input[pan_kmer_counts]} -o {params.output_dir} -c {threads} > {log} 2>&1
        """

rule normalize_kmer_counts:
    input:
        overall = "work/04.genotyping/overall_counts.tsv.gz",
        mapped = "work/04.genotyping/overall_mapped_kmers.tsv.gz"
    output:
        "work/04.genotyping/final.genotype.tsv",
        "work/04.genotyping/final.GAPIT.RData"
    params:
        script = SCRIPTS_DIR / "normalize_and_export.R",
        input_dir = "work/04.genotyping/",
        output_prefix = "work/04.genotyping/final"
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
        Rscript {params.script} --cm {input.mapped} --co {input.overall} -o {params.output_prefix} > {log} 2>&1
        """

rule genotyping_complete:
    input:
        "work/04.genotyping/final.genotype.tsv",
        "work/04.genotyping/final.GAPIT.RData"
    output:
        "work/04.genotyping/.done"
    log:
        "logs/04.genotyping/genotyping_complete.log"
    benchmark:
        "benchmarks/04.genotyping/genotyping_complete.tsv"
    shell:
        "touch {output}"
