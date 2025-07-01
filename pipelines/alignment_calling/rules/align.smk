

rule indexref:
    input:
        reference="results/ref/reference.fasta",
    output:
        indexedref="results/ref/reference.fasta.amb"
    log:
    message: "Indexing reference"
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        bwa index {input.reference}
        """


rule align_to_ref:
    input:
        R1 = 'results/rawdata/{sample}_R1.fastq.gz',
        R2 = 'results/rawdata/{sample}_R2.fastq.gz',
        indexedref="results/ref/reference.fasta.amb",
        ref="results/ref/reference.fasta"
    output:
        aligned = os.path.join(config['outdir'],'align/{sample}_aligned.bam') # .bam file output of aligned and sorted reads
    params:
        ref=config['reference']
    resources:
        runtime=600,
        mem_mb=lambda wc, input: max(2 * input.size_mb, 4000),
	    cores=4,
    log:
        stdout="logs/bwa/{sample}_aln.out",
        stderr="logs/bwa/{sample}_aln.err"
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    message: "Aligning reads to ref and sorting output for sample {wildcards.sample}"
    shell:
        """
        echo "Aligning reads for {wildcards.sample} to {input.ref}\n"
        echo 'bwa mem -o {output.aligned} {input.ref} {input.R1} {input.R2} \n' 
        bwa mem -t {resources.cores} {input.ref} {input.R1} {input.R2} \
            1> {output.aligned}  2> {log.stderr}
        """

rule primerclip:
    input:
        aligned = 'results/align/{sample}_aligned.bam',
        primers = config['primers']
    output:
        aln_trimmed='results/align/{sample}_aln_itrim.bam'
    group: "align"
    resources:
            mem_mb=8000,
            runtime=600,
    log:
        stdout="logs/ivar/{sample}_trim.out",
        stderr="logs/ivar/{sample}_trim.err"
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    message: "QC and soft-clipping primers using iVar"
    shell:
        """
        ivar trim -i {input.aligned} -b {input.primers} -p {output.aln_trimmed} -e > {log.stdout} 2> {log.stderr}
        """

rule sort:
    input:
        aln_trimmed= os.path.join(config['outdir'],'align/{sample}_aln_itrim.bam')
    output:
        aln_itrim_sorted=os.path.join(config['outdir'],'align/{sample}_aln_itrim_sorted.bam')
    resources:
        mem_mb=16000,
        runtime=180,
	cores=4,
    group: "align"
    log:
        stdout="logs/align/{sample}_sort.out",
        stderr="logs/align/{sample}_sort.err"
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    message: "Sorting and indexing reads"
    shell:
        """
        samtools sort {input.aln_trimmed} -@ {resources.cores} -o {output.aln_itrim_sorted} > {log.stdout} 2> {log.stderr}
        samtools index {output.aln_itrim_sorted} >> {log.stdout} 2>> {log.stderr}
        """

rule subsample:
    input:
        aln_itrim_sorted= os.path.join(config['outdir'],'align/{sample}_aln_itrim_sorted.bam')
    output:
        subsamp=os.path.join(config['outdir'],'align/{sample}_subsamp.bam')
    resources:
        mem_mb=8000,
        runtime=180,
    group: "align"
    params:
        subsample=config['subsample']
    log:
        stdout="logs/align/{sample}_sub.out",
        stderr="logs/align/{sample}_sub.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    message: "Selecting subsample of reads from {wildcards.sample}: {params.subsample}"
    shell:
        """
        samtools view -b -s {params.subsample} {input.aln_itrim_sorted} -o {output.subsamp} 1> {log.stdout} 2> {log.stderr}
        """
