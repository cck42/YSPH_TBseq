

rule bwa_index:
    input:
        fasta = "results/refs/{target}.fasta"
    output:
        bwt = "results/refs/{target}.fasta.bwt"
    log:
        log = "logs/bwa/index_{target}.log"
    params:
        idx_script = "scripts/indexer.sh",
    resources:
        mem_mb="40G",
        cpus_per_task=4,
        runtime=300
    container: "docker://sethnr/pgcoe_anypipe:0.01"
    shell:"""
        {params.idx_script} \
                -C {resources.cpus_per_task} \
                {input.fasta} >> {log.log} 2>&1
            """


rule bwa_align:
    input:
        R1 = 'results/rawdata/{sample}_R1.fastq.gz',
        R2 = 'results/rawdata/{sample}_R2.fastq.gz',
        indexedref="results/refs/{target}.fasta.bwt"
    output:
        aligned = temporary('results/align/{sample}_{target}_unsort.bam'),
        #aligned = 'results/align/{sample}_{target}_unsort.bam',
    group: "aligngroup"
    params:
        ref="results/refs/{target}.fasta",
    resources:
        runtime=180,
        mem_mb=lambda wc, input: max(2 * input.size_mb, 8000),
        cpus_per_task=4,
        cores=1,
    log:
        stderr="logs/align/bwa_mem_{sample}_{target}.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        echo "Aligning reads for {wildcards.sample} to {params.ref}\n"
        echo 'bwa mem -t {resources.cpus_per_task} {params.ref} {input.R1} {input.R2} \n  
                 -o {output.aligned}  \n' 
        bwa mem -t {resources.cpus_per_task} {params.ref} {input.R1} {input.R2} \
            -o {output.aligned}  2> {log.stderr}
        """

rule flagstat:
    input:
        aligned = 'results/align/{sample}_{target}_unsort.bam',
    output:
        flagstat = 'results/align/{sample}_{target}.flagstats'
    params:
        ref="results/refs/{target}.fasta",
        sleeplen=60,
    group: "aligngroup"
    resources:
        runtime=30,
        mem_mb=4000,
        cores=1,
    log:
        stderr="logs/align/flagstat_{sample}_{target}.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        samtools quickcheck {input.aligned} 2>&1 > {log.stderr}
        if [[ $? >0 ]]; then
            echo "sleeping {params.sleeplen}"
            sleep {params.sleeplen}
        else:
           echo "quickcheck good $?"
        fi

        samtools flagstat -@ {resources.cores} -O tsv   {input.aligned} > {output.flagstat}

        """

rule remove_unaligned:
    input:
        unsorted = 'results/align/{sample}_{target}_unsort.bam',
    output:
        aligned = temporary('results/align/{sample}_{target}_aligned.bam'),
    params:
        ref="results/refs/{target}.fasta",
        sleeplen=60,
    group: "aligngroup"
    resources:
        runtime=30,
        mem_mb=4000,
        cores=1,
    log:
        stderr="logs/align/flagstat_{sample}_{target}.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        samtools view -b -F 4 -F 2048 -o {output.aligned} {input.unsorted}
        """

rule sam_sort:
    input:
        aligned = 'results/align/{sample}_{target}_aligned.bam'
    output:
        sorted = 'results/align/{sample}_{target}.bam',
        idx = 'results/align/{sample}_{target}.bam.bai'
    resources:
        mem_mb=16000,
        runtime=40,
        cores=4,
    group: "aligngroup"
    log:
        stderr="logs/align/sort_{sample}_{target}.err"
    message: "Sorting and indexing reads"
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        samtools sort -@ {resources.cores} -o {output.sorted} {input.aligned} 2>&1 > {log.stderr}
        samtools index {output.sorted} 2>&1 >> {log.stderr}
        """




rule align_host:
    input:
        R1 = 'results/rawdata/{sample}_R1.fastq.gz',
        R2 = 'results/rawdata/{sample}_R2.fastq.gz',
        indexedref="results/refs/{host}.fasta.bwt"
    output:
        aligned=temporary('results/clean/{sample}_{host}_aligned.sam'),
        flagstat='results/clean/{sample}_{host}.flagstats',
    params:
        ref="results/refs/{host}.fasta",
        sleeplen=60,
    #group: "cleangroup"
    resources:
        runtime=600,
        mem_mb=lambda wc, input: max(2 * input.size_mb, 8000),
        cpus_per_task=4,
    log:
        stderr="logs/clean/{sample}_{host}_aligned.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        bwa mem -t {resources.cpus_per_task} -o {output.aligned} {params.ref} {input.R1} {input.R2} 2> {log.stderr}

        samtools flagstat -@ {resources.cpus_per_task} -O tsv  {output.aligned} > {output.flagstat}
        """


rule split_clean_dirty:
    input:
        aligned='results/clean/{sample}_{host}_aligned.sam',
        indexedref="results/refs/{host}.fasta.bwt"
    output:
        unaligned = 'results/clean/{sample}_{host}_cleaned.bam',
        flagstat='results/clean/{sample}_{host}_cleaned.flagstats',
        human='results/clean/{sample}_{host}.bam',
    params:
        ref="results/refs/{host}.fasta",
        mapqual=40,
    #group: "cleangroup"
    resources:
        runtime=180,
        mem_mb=lambda wc, input: max(2 * input.size_mb, 8000),
        cpus_per_task=1,
        cores=1,
    log:
        stderr="logs/clean/{sample}_{host}_cleaned.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        #samtools view -b -f 4 -o  {input.aligned} 
        samtools view -@ {resources.cpus_per_task} -b -F 4 -q {params.mapqual} \
             -U {output.unaligned} -o {output.human} {input.aligned}

        samtools flagstat -@ {resources.cpus_per_task} -O tsv {output.unaligned} > {output.flagstat}
        """

rule bam2fastq:
    input:
        unaligned = 'results/clean/{sample}_{host}_cleaned.bam',
    output:
        R1 = 'results/clean/{sample}_{host}_R1.fastq.gz',
        R2 = 'results/clean/{sample}_{host}_R2.fastq.gz',
    #group: "cleangroup"
    resources:
        runtime=120,
        mem_mb=lambda wc, input: max(2 * input.size_mb, 8000),
        cpus_per_task=1,
        cores=1,
    log:
        stderr="logs/clean/{sample}_{host}_cleaned.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        samtools fastq -1 {output.R1} -2 {output.R2}  {input.unaligned}
        """