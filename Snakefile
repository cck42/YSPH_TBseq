import pandas as pd
import os

configfile: "config/config.yaml"
include: "rules/common.py"

#RUNDATE=config['rundate']
SAMPLES = set(pd.read_table(config['samples'],header=0,dtype=str,delimiter=',').loc[lambda df: df['NGS_Date'] == '2024-02-12', 'Seq_ID'].tolist())
#SAMPLES = set(pd.read_table(config['samples'],header=0,dtype=str,delimiter=',')['Seq_ID']) #change from 'strain' as applicable
#SAMPLES=config['samples'] #un-comment for testing with single samples
#PATHS=set(pd.read_table(config['seq_data'],dtype=str))
#REFERENCE=config['reference']
# HOMEDIR=config['homedir'] #only if you want to set to something other than where this will be run from
#SCRATCH=config['scratchdir']
#PRIMERS=config['primers']
#OUTDIR=config['outdir']
READ_SUBSAMPLES=set(config['read_subsamples'])


rule all:
    input:
        coveragestats=expand(os.path.join(config['outdir'],'{sample}/{sample}_coveragestats_sub{subsample}.csv'),sample=SAMPLES,subsample=READ_SUBSAMPLES),
        vcf_filt=os.path.join(config['outdir'],config['runname'],"variants_filt.vcf.gz"),
        consensuses=expand(os.path.join(config['outdir'],'{sample}/{sample}_consensus.fa'),sample=SAMPLES),
        #snippyout=expand(os.path.join(config['outdir'],'{sample}/snippy_{sample}','snps.vcf'),sample=SAMPLES)
        combinedcoverage=os.path.join(config['outdir'],'TB001/combinedcoverage.tsv')
    
rule indexref:
    input:
        reference=config['reference']
    output:
        indexedref="data/reference.fasta.amb"
    log:
    message: "Indexing H37Rv reference"
    shell:
        """
        bwa index {input.reference}
        """

rule align_to_ref:
    input:
        read_location=os.path.join(config['outdir'],'{sample}','Unaligned'),
        indexedref="data/reference.fasta.amb"
    params:
        ref=config['reference']
    output:
        aligned = os.path.join(config['outdir'],'{sample}/{sample}_aligned.bam') # .bam file output of aligned and sorted reads
    log:
        stdout="logs/{sample}/aln.out",
        stderr="logs/{sample}/aln.err"
    message: "Aligning reads to H37Rv and sorting output for sample {wildcards.sample}"
    shell:
        """
        echo "Aligning reads for {wildcards.sample} to H37Rv reference strain\n"
        echo 'bwa mem {params.ref} {input.read_location}/*R1* {input.read_location}/*R2* | samtools view -b -F 4 -F 2048 | samtools sort -o {wildcards.sample}_aln.bam\n' 
        bwa mem {params.ref} {input.read_location}/*R1* {input.read_location}/*R2* | samtools view -b -F 4 -F 2048 | samtools sort -o {output.aligned} >> {log.stdout} 2>> {log.stderr}
        """

rule primerclip:
    input:
        aligned = os.path.join(config['outdir'],'{sample}/{sample}_aligned.bam'),
        primers = config['primers']
    output:
        aln_trimmed= os.path.join(config['outdir'],'{sample}/{sample}_aln_trimmed.bam')
    log:
        stdout="logs/{sample}/trim.out",
        stderr="logs/{sample}/trim.err"
    message: "QC and soft-clipping primers using iVar"
    shell:
        """
        ivar trim -i {input.aligned} -b {input.primers} -p {output.aln_trimmed} -e > {log.stdout} 2> {log.stderr}
        """

rule sort:
    input:
        aln_trimmed= os.path.join(config['outdir'],'{sample}/{sample}_aln_trimmed.bam')
    output:
        aln_trimmed_sorted=os.path.join(config['outdir'],'{sample}/{sample}_aln_trimmed_sorted.bam')
    log:
        stdout="logs/{sample}/sort.out",
        stderr="logs/{sample}/sort.err"
    message: "Sorting and indexing reads"
    shell:
        """
        samtools sort {input.aln_trimmed} -o {output.aln_trimmed_sorted} > {log.stdout} 2> {log.stderr}
        samtools index {output.aln_trimmed_sorted} >> {log.stdout} 2>> {log.stderr}
        """

rule coverageinfo:
    input:
        aln_trimmed_sorted= os.path.join(config['outdir'],'{sample}/{sample}_aln_trimmed_sorted.bam')
    output:
        coveragestats=os.path.join(config['outdir'],'{sample}/{sample}_coveragestats_sub{subsample}.csv'),
        subsamp=os.path.join(config['outdir'],'{sample}/{sample}_subsamp{subsample}.bam')
    params:
        maxdepth=0,
        minmapqual=60,
        minbasequal=13
    log:
        stdout="logs/{sample}/coverageinfo_sub{subsample}.out",
        stderr="logs/{sample}/coverageinfo_sub{subsample}.err",
    message: "Computing coverage of reference for sample {wildcards.sample} with a read subsample fraction of {wildcards.subsample}"
    shell:
        """
        samtools view -b -s {wildcards.subsample} {input.aln_trimmed_sorted} -o {output.subsamp} 1> {log.stdout} 2> {log.stderr}
        #echo -e "$(echo -e '{wildcards.sample}\t{wildcards.subsample}\t')\t" 1> {output.coveragestats} 2>> {log.stderr}
        #samtools coverage -d {params.maxdepth} -q {params.minmapqual} -Q {params.minbasequal} --no-header {output.subsamp} 1>> {output.coveragestats} 2>> {log.stderr}
        printf '\\n%s\\t%s\\t' {wildcards.sample} {wildcards.subsample} `samtools coverage -d {params.maxdepth} -q {params.minmapqual} -Q {params.minbasequal} --no-header {output.subsamp}` > {output.coveragestats} 2> {log.stderr}
        """

rule combinedcoverage:
    input:
        coveragestats=expand(os.path.join(config['outdir'],'{sample}/{sample}_coveragestats_sub{subsample}.csv'),sample=SAMPLES,subsample=READ_SUBSAMPLES)
    output:
        combinedcoverage=os.path.join(config['outdir'],'TB001/combinedcoverage.tsv')
    shell:
        """
        {{ echo -e "sample\tsubsample\trname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq";cat {input.coveragestats}; }} > {output.combinedcoverage}
        """

rule indexbam:
    input:
        bam = '{samplename}.bam'
    output:
        indexedbam = '{samplename}.bam.bai'
    log:
        stdout="logs/bamindex/{samplename}.out",
        stderr="logs/bamindex/{samplename}.err"
    message: "Indexing bam file: {wildcards.samplename}"
    shell:
        """
        samtools index {input.bam} > {log.stdout} 2> {log.stderr}
        """

rule variantcalling:
    input:
        tocall=expand(os.path.join(config['outdir'],'{sample}/{sample}_aln_trimmed_sorted.bam'),sample=SAMPLES),
        indexed=expand(os.path.join(config['outdir'],'{sample}/{sample}_aln_trimmed_sorted.bam.bai'),sample=SAMPLES)
    params:
        ref=config['reference']
    output:
        vcf=os.path.join(config['outdir'],config['runname'],"variants_unfilt.vcf.gz")
    log:
        stdout="logs/variantcalling.out",
        stderr="logs/variantcalling.err"
    message: "Calling variants for all samples"
    shell:
        """
        bcftools mpileup -Ou -o variants.bcf -f {params.ref} {input.tocall} 1> {log.stdout} 2> {log.stderr}
        bcftools call --ploidy 1 -vcO z -o {output.vcf} variants.bcf 1>> {log.stdout} 2>> {log.stderr}
        tabix -p vcf {output.vcf} 1>> {log.stdout} 2>> {log.stderr}
        """

rule variantfilter:
    input:
        vcf={rules.variantcalling.output.vcf}
    output:
        vcf_filt=os.path.join(config['outdir'],config['runname'],"variants_filt.vcf.gz")
    log:
        stdout="logs/variantfilter.out",
        stderr="logs/variantfilter.err"
    shell:
        """
        bcftools filter -O z -o {output.vcf_filt} -i 'QUAL>10 & DP>10' {input.vcf} 1> {log.stdout} 2> {log.stderr}
        """

rule consensus:
    input:
        tocall=os.path.join(config['outdir'],'{sample}/{sample}_aln_trimmed_sorted.bam'),
        indexed=os.path.join(config['outdir'],'{sample}/{sample}_aln_trimmed_sorted.bam.bai')
    params:
        ref=config['reference'],
        threshold=0.2,
        depth=20,
        prefix=os.path.join(config['outdir'],'{sample}/{sample}_consensus')
    output:
        consensus=os.path.join(config['outdir'],'{sample}/{sample}_consensus.fa')
    log:
        stdout="logs/{sample}/consensus.out",
        stderr="logs/{sample}/consensus.err"
    shell:
        """
        samtools mpileup -aa -A -d 0 -Q 0 -f {params.ref} {input.tocall} | ivar consensus -t {params.threshold} -m {params.depth} -p {params.prefix} -i {wildcards.sample} 1> {log.stdout} 2> {log.stderr}
        """

rule snippy:
    input:
        read_location=os.path.join(config['outdir'],'{sample}','Unaligned'),
        indexedref="data/reference.fasta.amb"
    params:
        ref=config['reference']
    output:
        snippyout = os.path.join(config['outdir'],'{sample}/snippy_{sample}','snps.vcf') # .bam file output of aligned and sorted reads
    log:
        stdout="logs/{sample}/snippy.out",
        stderr="logs/{sample}/snippy.err"
    message: "Running snippy on sample: {wildcards.sample}"
    shell:
        """
        snippy \
            --ref {params.ref} \
            --R1 {input.read_location}/*R1* \
            --R2 {input.read_location}/*R2* \
            --outdir results/{wildcards.sample}/snippy_{wildcards.sample} \
            --force 1> {log.stdout} 2> {log.stderr}
        """
