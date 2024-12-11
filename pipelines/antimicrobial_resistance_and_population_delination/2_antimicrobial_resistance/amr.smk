import pandas as pd 

configfile: "config/SP_read_aln.yaml"

samples_df = pd.read_csv("tsv/SP_czidseqs.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
CONTIGS = {row.sample_id: row.contig_path for row in samples_df.itertuples()}

rule all:
    input:
        expand("results/{genera}/abricate/{sample}/card_abricate.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/abricate/{sample}/megares_abricate.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/abricate/{sample}/argannot_abricate.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/abricate/{sample}/vfdb_abricate.tsv", sample=SAMPLES, genera=config["genera"]), 
        expand("results/{genera}/mlst/{sample}_mlst.csv", genera=config["genera"], sample=SAMPLES)

rule abricate:
    """
    Run in-silico AMR prediction
    """
    input:
        "assemblies/{sample}.fasta"
    output:
        out1="results/{genera}/abricate/{sample}/card_abricate.tsv",
        out2="results/{genera}/abricate/{sample}/megares_abricate.tsv",
        out3="results/{genera}/abricate/{sample}/argannot_abricate.tsv",
        out4="results/{genera}/abricate/{sample}/vfdb_abricate.tsv"
    params:
        genera=config["genera"]
    shell:
        """
        ml miniconda
        conda activate /home/flg9/.conda/envs/abricate
        abricate {input} --minid 75 --mincov 75 --db card > {output.out1}
        abricate {input} --minid 75 --mincov 75 --db megares > {output.out2}
        abricate {input} --minid 75 --mincov 75 --db argannot > {output.out3}
        abricate {input} --minid 75 --mincov 75 --db vfdb > {output.out4}
        """

rule mlst:
    """
    In-silico MLST assignment
    """
    input:
        "assemblies/{sample}.fasta"
    output:
        csv = "results/{genera}/mlst/{sample}_mlst.csv"
    shell:
        """
        ml miniconda 
        conda activate /home/flg9/.conda/envs/mlst
        mlst --csv {input} > {output.csv}
        """