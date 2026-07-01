def get_genome_file(wildcards):
    return genome_map[wildcards.genome_id]

rule complete_single_copy:
    input:
        "work/01.single_copy/single_copy.fa",
        "work/01.single_copy/single_copy.fa.fai"
    output:
        "work/01.single_copy/.done"
    shell:
        "touch {output}"

rule copy_genome:
    input:
        get_genome_file
    output:
        "work/01.single_copy/{genome_id}.fasta"
    shell:
        "cp {input} {output}"

rule genome_index:
    input:
        "work/01.single_copy/{genome_id}.fasta"
    output:
        "work/01.single_copy/{genome_id}.fasta.fai"
    conda:
        ENV_DIR / "single_copy.yaml"
    log:
        "logs/01.single_copy/{genome_id}.fasta.fai.log"
    shell:
        "samtools faidx {input}"

rule count_kmer:
    input:
        "work/01.single_copy/{genome_id}.fasta",
        "work/01.single_copy/{genome_id}.fasta.fai"
    output:
        temp("work/01.single_copy/{genome_id}.kmerCount.txt")
    conda:
        ENV_DIR / "single_copy.yaml"
    log:
        "logs/01.single_copy/{genome_id}.kmercount.log"
    resources:
        memory = config["max_mem"]
    params:
        mincount=config["min_kmer_count"],
        kmer_size=config["kmer_size"]
    shell:
        "kmercountexact.sh -Xmx{resources.memory}G mincount={params.mincount} in={input[0]} \
         out={output} k={params.kmer_size} 2> {log}"

rule mask_genome:
    input:
        "work/01.single_copy/{genome_id}.fasta",
        "work/01.single_copy/{genome_id}.kmerCount.txt"
    output:
        temp("work/01.single_copy/{genome_id}.masked.fa"),
        temp("work/01.single_copy/{genome_id}.masked.bbduk.out")
    params:
        kmer_size = config["kmer_size"]
    resources:
        memory = config["max_mem"]
    log:
        "logs/01.single_copy/{genome_id}.masked.bbduk.log"
    conda:
        ENV_DIR / "single_copy.yaml"
    shell:
        "bbduk.sh kmask='N' k={params.kmer_size} -Xmx{resources.memory}G in={input[0]} ref={input[1]} \
         out={output[0]} > {output[1]} 2> {log}"

rule index_masked_genome:
    input:
        "work/01.single_copy/{genome_id}.masked.fa"
    output:
        temp("work/01.single_copy/{genome_id}.masked.fa.fai")
    conda:
        ENV_DIR / "single_copy.yaml"
    log:
        "logs/01.single_copy/{genome_id}.masked.fa.fai.log"
    shell:
        "samtools faidx {input}"

rule extract_unmasked_regions:
    input:
        fasta="work/01.single_copy/{genome_id}.masked.fa",
        fai="work/01.single_copy/{genome_id}.masked.fa.fai"
    output:
        bed="work/01.single_copy/{genome_id}.noGaps.bed",
        bit=temp("work/01.single_copy/{genome_id}.masked.bit"),
        gaps_bed=temp("work/01.single_copy/{genome_id}.masked_gaps.bed"),
        genome_tmp=temp("work/01.single_copy/{genome_id}.masked.genome.tmp"),
    conda:
        ENV_DIR / "single_copy.yaml"
    log:
        "logs/01.single_copy/{genome_id}.unmasked.log"
    shell:
        r"""
            (
        # Create FASTA index file
        samtools faidx {input.fasta}

        # Convert fasta to twoBit
        faToTwoBit -long {input.fasta} {output.bit}

        # Get gaps in BED format
        twoBitInfo -nBed {output.bit} {output.gaps_bed}

        # Prepare genome file
        cut -f1,2 {input.fai} > {output.genome_tmp}

        # Find non-gap regions as final BED output
        bedtools complement -g {output.genome_tmp} -i {output.gaps_bed} > {output.bed}
        ) > {log} 2>&1
        """

rule extract_single_copy_regions:
    input:
        bed="work/01.single_copy/{genome_id}.noGaps.bed",
        fa="work/01.single_copy/{genome_id}.fasta"
    output:
        fa="work/01.single_copy/{genome_id}.single_copy.fa",
        bed="work/01.single_copy/{genome_id}.single_copy.bed"
    params:
        minlen=config["min_length"],
        mergegap=config["merge_gap"],
    conda:
        ENV_DIR / "single_copy.yaml"
    log:
        "logs/01.single_copy/{genome_id}.single_copy.log"
    shell:
        r"""
            (
            awk -v minlen={params.minlen} '$3-$2 >= minlen' {input.bed} | bedtools merge -d {params.mergegap} > {output.bed}
            awk -v name={wildcards.genome_id} '{{print $0"\t"name"|"$1":"$2+1"-"$3}}' {output.bed} \
                | bedtools getfasta -bed - -fi {input.fa} -fo /dev/stdout -name \
                | sed 's/::.*$//' > {output.fa}
        ) > {log} 2>&1
        """

rule concat_single_copy_regions:
    input:
        expand("work/01.single_copy/{genome_id}.single_copy.fa", genome_id=genome_df.index)
    output:
        "work/01.single_copy/single_copy.fa"
    shell:
        r"""
        cat {input} > {output}
        """

rule index_single_copy_fasta:
    input:
        "work/01.single_copy/single_copy.fa"
    output:
        "work/01.single_copy/single_copy.fa.fai"
    conda:
        ENV_DIR / "single_copy.yaml"
    shell:
        r"""
        samtools faidx {input}
        """
