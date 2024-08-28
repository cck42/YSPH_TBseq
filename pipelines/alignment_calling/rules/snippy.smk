rule snippy:
    input:
        #read_location=os.path.join(config['outdir'],'{sample}','Unaligned'),
        #indexedref="data/reference.fasta.amb"
        read_location=os.path.join(config['readdir'],'{sample}'),
        indexedref="data/reference.fasta.amb"
    params:
        ref=config['reference']
    container: "docker://staphb/snippy"
    resources:
            mem_mb=32000,
            runtime=240,
	    cores=4,
    output:
        vcf = os.path.join(config['outdir'],'snippy/snippy_{sample}','snps.vcf'),
        tab = os.path.join(config['outdir'],'snippy/snippy_{sample}','snps.tab')
    log:
        stdout="logs/snippy/{sample}.out",
        stderr="logs/snippy/{sample}.err"
    message: "Running snippy on sample: {wildcards.sample}"
    shell:
        """
        snippy \
            --ref {params.ref} \
            --cpus {resources.cores} \
            --R1 {input.read_location}/*R1* \
            --R2 {input.read_location}/*R2* \
            --outdir results/snippy/snippy_{wildcards.sample} \
            --force 1> {log.stdout} 2> {log.stderr}
        """

rule snippy_core:
    input:
        snippyout=expand(config['outdir']+'/snippy/snippy_{sample}/snps.tab',sample=SAMPLES),
        indexedref="data/reference.fasta.amb"
    params:
        ref=config['reference'],
        prefix=config['outdir']+'/snippy/core',
        indirs=expand(config['outdir']+'/snippy/snippy_{sample}',sample=SAMPLES)
    container: "docker://staphb/snippy"
    resources:
            mem_mb=32000,
            runtime=240,
    output:
        coreout = expand(config['outdir']+'/snippy/core.{suffix}',suffix=['aln','tab','txt','vcf']) 
    log:
        stdout="logs/"+config['runname']+"/core.out",
        stderr="logs/"+config['runname']+"/core.err"
    message: "Running snippy core"
    shell:
        """
        snippy-core --ref {params.ref}  --prefix {params.prefix} \
            {params.indirs} \
            1> {log.stdout} 2> {log.stderr}
        """
