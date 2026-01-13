rule cp_local_fq:
    input:
        read_location=os.path.join(config['readdir'],'{sample}'),
    output:
        R1 = 'results/rawdata/{sample}_R1.fastq.gz',
        R2 = 'results/rawdata/{sample}_R2.fastq.gz'
    log:
        log = "logs/datacopy/{sample}.log"
    group: "cplocal"
    resources:
        mem_mb="2G",
        cpus_per_task=1,
        runtime=180
    container: None
    shell:"""
        #find {input.read_location}

        R1=` find {input.read_location} | grep _R1_ | {{ grep .gz$ || true;}}`
        R2=` find {input.read_location} | grep _R2_ | {{ grep .gz$ || true;}}`

        #check for quip and convert to gz if needed
        ml Quip 2>>{log.log}
            QP1=` find {input.read_location} | grep _R1_ | {{ grep .qp$ || true;}}`
            QP2=` find {input.read_location} | grep _R2_ | {{ grep .qp$ || true;}}`

            # If R1 gz file is missing but R1 qp file exists, convert both R1 and R2 qp to gz in place
            if [ -z "$R1" ] && [ -f "$QP1" ]; then
                quip -dc "$QP1" | gzip > "${{QP1/.qp/.gz}}"
                quip -dc "$QP2" | gzip > "${{QP2/.qp/.gz}}"
                R1=` find {input.read_location} | grep _R1_ | grep .gz$`
                R2=` find {input.read_location} | grep _R2_ | grep .gz$`
            fi

        echo cp $R1 {output.R1}
        cp $R1 {output.R1}
        
        echo cp $R2 {output.R2}
        cp $R2 {output.R2}
        
        """


rule cp_local_data:
    input:
        ref=config['reference'],
        gff=config['gff'],
        primers=config['primers'],
        amplicons=config['amplicons']
    output:
        ref = 'results/ref/reference.fasta',
        gff = 'results/ref/genes.gff',
        amplicons = 'results/ref/amplicons.bed',
        primers = 'results/ref/primers.bed',
    log:
        log = "logs/datacopy/copyrefs.log"
    group: "cplocal"
    resources:
        mem_mb="2G",
        cpus_per_task=1,
        runtime=60
    params:
        targets = TARGETS
    container: None
    shell:"""
        cp {input.ref} {output.ref}
        cp {input.gff} {output.gff}
        cp {input.primers} {output.primers}
        cp {input.amplicons} {output.amplicons}
        """

