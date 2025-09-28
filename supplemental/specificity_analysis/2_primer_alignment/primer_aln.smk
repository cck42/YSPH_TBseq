configfile: "config/TB_primer_aln_reps.yaml"

rule all:
    input:
        expand("results/{genera}/specificity_analysis/primer_aln/split_primers/fwd_primers.bed", genera=config["genera"]), 
        expand("results/{genera}/specificity_analysis/primer_aln/split_primers/rev_primers.bed", genera=config["genera"]), 
        expand("results/{genera}/specificity_analysis/primer_aln/get_fasta/fwd.scheme.primer.fasta", genera=config["genera"]),
        expand("results/{genera}/specificity_analysis/primer_aln/get_fasta/rev.scheme.primer.fasta", genera=config["genera"]),
        expand("results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.1.bt2",sample=config["samples"], genera=config["genera"]), 
        expand("results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.2.bt2",sample=config["samples"], genera=config["genera"]), 
        expand("results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.3.bt2",sample=config["samples"], genera=config["genera"]), 
        expand("results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.4.bt2",sample=config["samples"], genera=config["genera"]), 
        expand("results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.rev.1.bt2",sample=config["samples"], genera=config["genera"]), 
        expand("results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.rev.2.bt2",sample=config["samples"], genera=config["genera"]), 
        expand("results/{genera}/specificity_analysis/align_primers/{sample}_fwd_aln.bam", sample=config["samples"], genera=config["genera"]),
        expand("results/{genera}/specificity_analysis/align_primers/{sample}_rev_aln.bam", sample=config["samples"], genera=config["genera"]),
        expand("results/{genera}/specificity_analysis/aln_loc/{sample}_fwd.depth", sample=config["samples"], genera=config["genera"]), 
        expand("results/{genera}/specificity_analysis/aln_loc/{sample}_rev.depth", sample=config["samples"], genera=config["genera"])

rule split_primers:
    """
    Use awk to extract fwd and rev primers from BED file
    """
    input:
        bed=config["bed"]
    output:
        fwd_pr="results/{genera}/specificity_analysis/primer_aln/split_primers/fwd_primers.bed", 
        rev_pr="results/{genera}/specificity_analysis/primer_aln/split_primers/rev_primers.bed"
    shell:
        """
        awk '$6 == "+" {{print}}' {input} > {output.fwd_pr}
        awk '$6 == "-" {{print}}' {input} > {output.rev_pr}
        """

rule get_fasta:
    """
    Convert BED files to FASTA files for downstream alignment 
    """
    input:
        fwd_pr="results/{genera}/specificity_analysis/primer_aln/split_primers/fwd_primers.bed", 
        rev_pr="results/{genera}/specificity_analysis/primer_aln/split_primers/rev_primers.bed"
    output:
        bedfasta_fwd = "results/{genera}/specificity_analysis/primer_aln/get_fasta/fwd.scheme.primer.fasta", 
        bedfasta_rev = "results/{genera}/specificity_analysis/primer_aln/get_fasta/rev.scheme.primer.fasta"
    log:
        stdout = "logs/{genera}/specificity_analysis/primer_aln/get_fasta/get_fasta.out",
        stderr = "logs/{genera}/specificity_analysis/primer_aln/get_fasta/get_fasta.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/biopython

        python \
        /vast/palmer/scratch/turner/flg9/GLab/specificity_analysis/2_primer_aln/scripts/get_indiv_fasta.py {input.fwd_pr} {input.rev_pr} {output.bedfasta_fwd} {output.bedfasta_rev} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule index_seqs:
    """
    Create index of rep sequences
    """
    input:
        queryseqs = lambda wildcards: config['samples'][wildcards.sample]
    output:
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.1.bt2",
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.2.bt2",
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.3.bt2",
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.4.bt2",
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.rev.1.bt2",
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.rev.2.bt2"
    params:
        genera=config["genera"],
        outdir="results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}"
    log:
        stdout = "logs/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/index_seqs.out",
        stderr = "logs/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/index_seqs.err"
    shell:
        """
        module unload miniconda
        module load Bowtie2/2.5.1-GCC-12.2.0

        bowtie2-build \
        -f {input.queryseqs} {params.outdir}/{wildcards.sample}_indexed_ref \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule align_primers:
    """
    Align primers to rep sequences
    """
    input:
        fwd = rules.get_fasta.output.bedfasta_fwd, 
        rev = rules.get_fasta.output.bedfasta_rev,
        idx = [
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.1.bt2",
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.2.bt2",
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.3.bt2",
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.4.bt2",
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.rev.1.bt2",
        "results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}/{sample}_indexed_ref.rev.2.bt2" 
        ]
    output:
        fwd_aln = "results/{genera}/specificity_analysis/align_primers/{sample}_fwd_aln.bam",
        rev_aln = "results/{genera}/specificity_analysis/align_primers/{sample}_rev_aln.bam"
    params:
        genera=config["genera"],
        dir="results/{genera}/specificity_analysis/primer_aln/index_seqs/{sample}"
    log:
        stdout = "logs/{genera}/specificity_analysis/align_primers/{sample}_align_primers.out",
        stderr = "logs/{genera}/specificity_analysis/align_primers/{sample}_align_primers.err"
    shell:
        """
        module unload miniconda
        module load Bowtie2/2.5.1-GCC-12.2.0
        module load SAMtools/1.21-GCC-12.2.0

        bowtie2 \
        -x {params.dir}/{wildcards.sample}_indexed_ref -f {input.fwd} -N 1 --mp 6,2 --rdg 5,1 --rfg 5,1 -L 15 -i S,1,0.50 \
        | samtools view -b -F 4 -F 2048 | samtools sort -o {output.fwd_aln} \
        1>> {log.stdout} 2>> {log.stderr}
        
        bowtie2 \
        -x {params.dir}/{wildcards.sample}_indexed_ref -f {input.rev} -N 1 --mp 6,2 --rdg 5,1 --rfg 5,1 -L 15 -i S,1,0.50 \
        | samtools view -b -F 4 -F 2048 | samtools sort -o {output.rev_aln} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule aln_loc:
    """
    Find where our primers are aligning by using samtools depth to give us a "depth" of reads
    that align to coordinates across our query genomes
    """
    input:
        fwd_primer_alns = rules.align_primers.output.fwd_aln,
        rev_primer_alns = rules.align_primers.output.rev_aln 
    output:
        fwd_primer_depth = "results/{genera}/specificity_analysis/aln_loc/{sample}_fwd.depth",
        rev_primer_depth = "results/{genera}/specificity_analysis/aln_loc/{sample}_rev.depth"
    log:
        stdout = "logs/{genera}/specificity_analysis/aln_loc/{sample}_aln_loc.out",
        stderr = "logs/{genera}/specificity_analysis/aln_loc/{sample}_aln_loc.err"
    shell:
        """
        module unload miniconda
        module load SAMtools/1.21-GCC-12.2.0

        samtools depth \
        -a {input.fwd_primer_alns} > {output.fwd_primer_depth}

        samtools depth \
        -a {input.rev_primer_alns} > {output.rev_primer_depth}
        """