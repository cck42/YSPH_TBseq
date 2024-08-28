configfile: "config/TB_primer_aln_reps.yaml"

rule all:
    input:
        expand("results/{genera}/primer_aln/fwd_primers.bed", genera=config["genera"]), 
        expand("results/{genera}/primer_aln/rev_primers.bed", genera=config["genera"]), 
        expand("results/{genera}/primer_aln/fwd.scheme.primer.fasta", genera=config["genera"]),
        expand("results/{genera}/primer_aln/rev.scheme.primer.fasta", genera=config["genera"]),
        expand("results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.1.bt2", sample=config["samples"], genera=config["genera"]),
        expand("results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.2.bt2", sample=config["samples"], genera=config["genera"]),
        expand("results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.3.bt2", sample=config["samples"], genera=config["genera"]),
        expand("results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.4.bt2", sample=config["samples"], genera=config["genera"]),
        expand("results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.1.bt2", sample=config["samples"], genera=config["genera"]),
        expand("results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.2.bt2", sample=config["samples"], genera=config["genera"]),
        expand("results/{genera}/primer_aln/bowtie_align/{sample}_fwd_aln.bam", sample=config["samples"], genera=config["genera"]),
        expand("results/{genera}/primer_aln/bowtie_align/{sample}_rev_aln.bam", sample=config["samples"], genera=config["genera"]),
        expand("results/{genera}/primer_aln/samtools_depth_indiv_primers/{sample}_fwd.depth", sample=config["samples"], genera=config["genera"]), 
        expand("results/{genera}/primer_aln/samtools_depth_indiv_primers/{sample}_rev.depth", sample=config["samples"], genera=config["genera"])

rule split_bed:
    """
    Use awk to extract fwd and rev primers 
    """
    input:
        bed=config["bed"]
    output:
        fwd_pr="results/{genera}/primer_aln/fwd_primers.bed", 
        rev_pr="results/{genera}/primer_aln/rev_primers.bed"
    shell:
        """
        awk '$6 == "+" {{print}}' {input} > {output.fwd_pr} 
        awk '$6 == "-" {{print}}' {input} > {output.rev_pr} 
        """

rule get_fasta:
    """
    Convert BED to FASTA for downstream alignment 
    """
    input:
        fwd_bed = "results/{genera}/primer_aln/fwd_primers.bed", 
        rev_bed = "results/{genera}/primer_aln/rev_primers.bed"
    output:
        bedfasta_fwd = "results/{genera}/primer_aln/fwd.scheme.primer.fasta", 
        bedfasta_rev = "results/{genera}/primer_aln/rev.scheme.primer.fasta"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python /vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/primer_alignment/scripts/get_indiv_fasta.py {input.fwd_bed} {input.rev_bed} {output.bedfasta_fwd} {output.bedfasta_rev}
        """

rule bowtie_build:
    """
    Create index of rep sequences
    """
    input:
        queryseqs = lambda wildcards: config['samples'][wildcards.sample]
    output:
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.1.bt2",
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.2.bt2",
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.3.bt2",
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.4.bt2",
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.1.bt2",
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.2.bt2"
    conda:
        "envs/read_aln.yaml"
    params:
        genera=config["genera"]
    log:
        "results/{genera}/primer_aln/logs/bowtie_build/{sample}_bowtie_build.log"
    shell:
        """
        bowtie2-build -f {input.queryseqs} results/{params.genera}/primer_aln/bowtie_index/{wildcards.sample}/{wildcards.sample}_indexed_ref > {log}
        """

rule bowtie_align:
    """
    Align primers to rep sequences
    """
    input:
        fwd = rules.get_fasta.output.bedfasta_fwd, 
        rev = rules.get_fasta.output.bedfasta_rev,
        idx = [
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.1.bt2",
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.2.bt2",
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.3.bt2",
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.4.bt2",
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.1.bt2",
        "results/{genera}/primer_aln/bowtie_index/{sample}/{sample}_indexed_ref.rev.2.bt2" 
        ]
    output:
        fwd_aln = "results/{genera}/primer_aln/bowtie_align/{sample}_fwd_aln.bam",
        rev_aln = "results/{genera}/primer_aln/bowtie_align/{sample}_rev_aln.bam"
    conda:
        "envs/read_aln.yaml"
    params:
        genera=config["genera"]
    log:
        "results/{genera}/primer_aln/logs/bowtie_align/bowtie_{sample}.log"
    shell:
        """
        bowtie2 -x results/{params.genera}/primer_aln/bowtie_index/{wildcards.sample}/{wildcards.sample}_indexed_ref -f {input.fwd} | samtools view -b -F 4 -F 2048 | samtools sort -o {output.fwd_aln} 2> {log}
        bowtie2 -x results/{params.genera}/primer_aln/bowtie_index/{wildcards.sample}/{wildcards.sample}_indexed_ref -f {input.rev} | samtools view -b -F 4 -F 2048 | samtools sort -o {output.rev_aln} 2> {log}
        """

rule fwd_rev_depth:
    """
    Calculate depth for fwd and rev primer alns separately
    """
    input:
        fwd_primer_alns = rules.bowtie_align.output.fwd_aln,
        rev_primer_alns = rules.bowtie_align.output.rev_aln 
    output:
        fwd_primer_depth = "results/{genera}/primer_aln/samtools_depth_indiv_primers/{sample}_fwd.depth",
        rev_primer_depth = "results/{genera}/primer_aln/samtools_depth_indiv_primers/{sample}_rev.depth"
    conda:
        "envs/read_aln.yaml"
    shell:
        """
        samtools depth -a {input.fwd_primer_alns} > {output.fwd_primer_depth}
        samtools depth -a {input.rev_primer_alns} > {output.rev_primer_depth}
        """