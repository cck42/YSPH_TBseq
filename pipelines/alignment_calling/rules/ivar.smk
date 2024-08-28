
rule pileup:
    input:
        tocall=os.path.join(config['outdir'],'align/{sample}_aln_itrim_sorted.bam'),
        indexed=os.path.join(config['outdir'],'align/{sample}_aln_itrim_sorted.bam.bai')
    params:
        ref=config['reference'],
        threshold=0.2,
        depth=20,
        prefix=os.path.join(config['outdir'],'ivar/{sample}_consensus')
    output:
        pileup=os.path.join(config['outdir'],'pileup/{sample}.pileup')
    log:
        stderr="logs/pileup/{sample}.err"
    shell:
        """
        samtools mpileup -aa -A -d 0 -Q 0 -f {params.ref} {input.tocall} 1> {output.pileup} 2> {log.stderr}
        """


rule consensus:
    input:
        pileup=os.path.join(config['outdir'],'pileup/{sample}.pileup')
    params:
        ref=config['reference'],
        threshold=0.2,
        depth=20,
        prefix=os.path.join(config['outdir'],'ivar/{sample}_consensus')
    output:
        consensus=os.path.join(config['outdir'],'ivar/{sample}_consensus.fa')
    log:
        stdout="logs/{sample}/consensus.out",
        stderr="logs/{sample}/consensus.err"
    shell:
        """
        cat {input.pileup} | ivar consensus -t {params.threshold} \
           -m {params.depth} -p {params.prefix} -i {wildcards.sample} \
           1> {log.stdout} 2> {log.stderr}
        """

rule ivariants:
    input:
        pileup=os.path.join(config['outdir'],'pileup/{sample}.pileup')
    params:
        ref=config['reference'],
        threshold=0.2,
        depth=20,
        qual=20,
        prefix=os.path.join(config['outdir'],'ivar/{sample}_ivariants')
    output:
        ivariants=os.path.join(config['outdir'],'ivar/{sample}_ivariants.tsv'),
    resources:
        mem_mb=8000,
        runtime=240,
    log:
        stderr="logs/ivar/{sample}_consensus.err"
    shell:
        """
        cat {input.pileup} | \
         ivar variants -q {params.qual} -r {params.ref} -t {params.threshold} \
           -m {params.depth} -p {params.prefix}  2>&1 >  {log.stderr}        
        """
