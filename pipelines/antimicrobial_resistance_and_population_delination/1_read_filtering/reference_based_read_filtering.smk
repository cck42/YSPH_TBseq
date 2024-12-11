import pandas as pd

configfile: "config/SP_read_aln.yaml"

samples_df = pd.read_csv("tsv/SP_seqs.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

rule all:
    input:
        expand("results/{genera}/ref_index/indexed_ref.0123", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.amb", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.ann", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.bwt.2bit.64", genera=config["genera"]),
        expand("results/{genera}/ref_index/indexed_ref.pac", genera=config["genera"]),
        expand("results/{genera}/alignments/{sample}/{sample}_aligned_sorted.bam", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/sorted_bams/{sample}/{sample}_qsort.bam", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/fastqs/{sample}/{sample}_R1.fq", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/fastqs/{sample}/{sample}_R2.fq", sample=SAMPLES, genera=config["genera"])

rule bwa_build:
    """
    Create index of OXC141
    """
    input:
        ref=config["ref"]
    output:
        "results/{genera}/ref_index/indexed_ref.0123",
        "results/{genera}/ref_index/indexed_ref.amb",
        "results/{genera}/ref_index/indexed_ref.ann",
        "results/{genera}/ref_index/indexed_ref.bwt.2bit.64",
        "results/{genera}/ref_index/indexed_ref.pac"
    params:
        genera=config["genera"]
    log:
        "results/{genera}/logs/bwa_build/build.log"
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        mkdir -p results/{params.genera}/ref_index/
        bwa-mem2 index {input.ref} -p results/{params.genera}/ref_index/indexed_ref > {log} 2>&1
        """

rule bwa:
    """
    Align raw reads to OXC141
    """
    input:
        fwd=lambda wildcards: READS[wildcards.sample]["r1"],
        rev=lambda wildcards: READS[wildcards.sample]["r2"],
        idx=[
            "results/{genera}/ref_index/indexed_ref.0123",
            "results/{genera}/ref_index/indexed_ref.amb",
            "results/{genera}/ref_index/indexed_ref.ann",
            "results/{genera}/ref_index/indexed_ref.bwt.2bit.64",
            "results/{genera}/ref_index/indexed_ref.pac"
        ]
    output:
        bam="results/{genera}/alignments/{sample}/{sample}_aligned_sorted.bam"
    params:
        genera=config["genera"]
    conda:
        "envs/read_aln.yaml"
    shell: 
        """
        bwa-mem2 mem -t 4 results/{params.genera}/ref_index/indexed_ref {input.fwd} {input.rev} > results/{params.genera}/alignments/{wildcards.sample}/{wildcards.sample}_aligned.sam
        samtools view -h -b -F 4 -F 2048 results/{params.genera}/alignments/{wildcards.sample}/{wildcards.sample}_aligned.sam | samtools sort -o {output.bam}
        """

rule sort_bams:
    """
    Sort by read names
    """
    input:
        "results/{genera}/alignments/{sample}/{sample}_aligned_sorted.bam"
    output:
        "results/{genera}/sorted_bams/{sample}/{sample}_qsort.bam"
    conda:
        "envs/read_aln.yaml"
    params:
        genera=config["genera"]
    shell:
        """
        samtools sort -n -o {output} {input}
        """

rule bamtofastq:
    """
    Filter out reads that aligned to OXC141 and convert to fastq
    """
    input:
        "results/{genera}/sorted_bams/{sample}/{sample}_qsort.bam"
    output:
        r1 = "results/{genera}/fastqs/{sample}/{sample}_R1.fq",
        r2 = "results/{genera}/fastqs/{sample}/{sample}_R2.fq"
    params:
        genera=config["genera"]
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        bedtools bamtofastq -i {input} -fq {output.r1} -fq2 {output.r2}
        """