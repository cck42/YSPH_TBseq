

# rule poppunkFile:
#     input:
#         R1s=expand('results/rawdata/{sample}_R1.fastq.gz',sample=SAMPLES),
#         R2s=expand('results/rawdata/{sample}_R2.fastq.gz',sample=SAMPLES)
#     output:
#         qfile="results/poppunk/query_file.txt",
#     params:
#         samples=SAMPLES,
#     run:
#         with open(output.qfile, 'w') as f:
#             for S in SAMPLES:
#                 print("{}\tresults/rawdata/{}_R1.fastq.gz\tresults/rawdata/{}_R2.fastq.gz".format(S,S,S),file=f)
            
# rule poppunk:
#     """
#     GPSC assignment
#     """
#     input:
#         queryseqs = "results/poppunk/query_file.txt",
#     output:
#         clusters="results/poppunk/GPSC_assignments_clusters.csv",
#         exclusters="results/poppunk/GPSC_assignments_external_clusters.csv",
#         fin="results/poppunk/complete"
#     resources:
#         mem_mb=64000,
#         runtime=240,
#         disk_mb=10000,
#         cores=12
#     container:"docker://staphb/poppunk"
#     params:
#         db=config['poppunkdbdir']
#     shell:
#         """
#         poppunk_assign --db {params.db} --serial \
#         --output results/poppunk --query {input.queryseqs} --threads 12

#         touch {output.fin}
#         """

rule punkone:
    """
    GPSC assignment
    """
    input:
        R1='results/rawdata/{sample}_R1.fastq.gz',
        R2='results/rawdata/{sample}_R2.fastq.gz'
    output:
        qfile="results/poppunk/{sample}/qfile.txt",
        clusters="results/poppunk/{sample}/{sample}_clusters.csv",
    resources:
        disk_mb=10000,
        mem_mb=16000,
        runtime=60,
        cores=1
    container:"docker://staphb/poppunk"
    params:
        db=config['poppunkdbdir']
    shell:
        """
        echo "{wildcards.sample}\t{input.R1}\t{input.R2}" > {output.qfile}

        poppunk_assign --db {params.db} \
        --output results/poppunk/{wildcards.sample} --query {output.qfile} --threads 1

        
        """

rule catpunk:
    input:
        clusters=expand('results/poppunk/{sample}/{sample}_clusters.csv',sample=SAMPLES),
    output:
        clusters="results/poppunk/clusters.csv",
        fin="results/poppunk/complete"
    params:
        samples=SAMPLES,
    shell:
        """
        cat {input.clusters} > {output.clusters}
        touch {output.fin}
        
        """

rule srstone:
    """
    GPSC assignment
    """
    input:
        R1='results/rawdata/{sample}_R1.fastq.gz',
        R2='results/rawdata/{sample}_R2.fastq.gz'
    output:
        clusters="results/srst/{sample}__mlst__{pathogen}__results.txt",
    resources:
        disk_mb=10000,
        mem_mb=16000,
        runtime=120,
        cores=1
    container:"docker://staphb/srst2"
    params:
        db=config['srstdb'],
        prefix="results/srst/{sample}",
    shell:
        """
        srst2 --output {params.prefix} \
          --input_pe {input.R1} {input.R2} \
          --forward _R1 --reverse _R2 \
          --mlst_db {params.db} \
          --mlst_definitions results/srstDB/profiles_csv \
          --mlst_delimiter '_'
        
        """

rule catsrst:
    input:
        clusters=expand("results/srst/{sample}__mlst__{pathogen}__results.txt",sample=SAMPLES,pathogen=config['pathogen']),
    output:
        clusters="results/srst/clusters.csv",
        fin="results/srst/complete"
    params:
        samples=SAMPLES,
    shell:
        """
        cat {input.clusters} > {output.clusters}
        touch {output.fin}
        
        """

# rule poppunkDB:
#     input:
#         queryseqs = "txt/qfile_contigs.txt"
#     output:
#         "results/poppunkDB/",
#     params:
#         db=config['poppunkdb']
#     shell:
#         """
#         mkdir results/poppunkDB/
#         cd results/poppunkDB/
#         wget {params.db}
#         """

