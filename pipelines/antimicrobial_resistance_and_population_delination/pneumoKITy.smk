import pandas as pd

configfile: "config/SP_read_aln.yaml"

samples_df = pd.read_csv("tsv/SP_seqs.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

rule all:
    input:
        expand("results/{genera}/pneumoKITy_pure/{sample}/pneumo_capsular_typing/{sample}_alldata.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/pneumoKITy_pure/{sample}/collated/Collated_result_data.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/pneumoKITy_mixed/{sample}/pneumo_capsular_typing/{sample}_alldata.csv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/pneumoKITy_mixed/{sample}/collated/Collated_result_data.csv", sample=SAMPLES, genera=config["genera"])

rule pneumoKITY_pure:
    """
    Run pneumoKITY pure on raw read SPn samples
    """
    input:
        fwd="results/{genera}/fastqs/{sample}/{sample}_R1.fq",
        rev="results/{genera}/fastqs/{sample}/{sample}_R2.fq"
    output:
        "results/{genera}/pneumoKITy_pure/{sample}/pneumo_capsular_typing/{sample}_alldata.csv",
        "results/{genera}/pneumoKITy_pure/{sample}/collated/Collated_result_data.csv"
    params:
        genera=config["genera"],
        outdir_pure="results/{genera}/pneumoKITy_pure/{sample}",
        collated_pure="results/{genera}/pneumoKITy_pure/{sample}/collated",
        threads=4
    shell:
        """
        ml miniconda
        conda activate /home/flg9/.conda/envs/pneumokity
        export PYTHONPATH=/home/flg9/.conda/envs/pneumokity/lib/python3.12/site-packages
        mkdir -p {params.collated_pure}
        python PneumoKITy/pneumokity.py pure -f {input.fwd} {input.rev} -t {params.threads} -o {params.outdir_pure} -c {params.collated_pure} -s {wildcards.sample} -n 10 -p 75
        """

rule pneumoKITY_mixed:
    """
    Run pneumoKITY mix on raw read SPn samples
    """
    input:
        fwd="results/{genera}/fastqs/{sample}/{sample}_R1.fq",
        rev="results/{genera}/fastqs/{sample}/{sample}_R2.fq"
    output:
        "results/{genera}/pneumoKITy_mixed/{sample}/pneumo_capsular_typing/{sample}_alldata.csv",
        "results/{genera}/pneumoKITy_mixed/{sample}/collated/Collated_result_data.csv"
    params:
        genera=config["genera"],
        outdir_mixed="results/{genera}/pneumoKITy_mixed/{sample}",
        collated_mixed="results/{genera}/pneumoKITy_mixed/{sample}/collated",
        threads=4
    shell:
        """
        ml miniconda
        conda activate /home/flg9/.conda/envs/pneumokity
        export PYTHONPATH=/home/flg9/.conda/envs/pneumokity/lib/python3.12/site-packages
        mkdir -p {params.collated_mixed}
        python PneumoKITy/pneumokity.py mix -f {input.fwd} {input.rev} -t {params.threads} -o {params.outdir_mixed} -c {params.collated_mixed} -s {wildcards.sample} -n 4 -p 75
        """