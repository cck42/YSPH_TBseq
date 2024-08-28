import pandas as pd

configfile: "config/TB_pangenome.yaml"

samples_df = pd.read_csv("tsv/TB_reps.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
SEQ = {row.sample_id: {"sample_id": row.sample_id, "seq": row.seq_path} for row in samples_df.itertuples()} 

rule all:
    input:
        expand("results/{genera}/prokka/{sample}/{sample}.gff", genera=config["genera"],sample=SAMPLES),
        expand("results/{genera}/roary/core_gene_alignment.aln", genera=config["genera"]),
        expand("results/{genera}/roary/gene_presence_absence.csv", genera=config["genera"]),
        expand("results/{genera}/roary/tree.newick", genera=config["genera"]), 
        expand("results/{genera}/fastani/fastani.out", genera=config["genera"])

rule prokka:
    """
    Obtain gff files needed for Roary
    """
    input:
        inputseqs=lambda wildcards: SEQ[wildcards.sample]["seq"]
    output:
        "results/{genera}/prokka/{sample}/{sample}.gff"
    params:
        target=config["target"], 
        genera=config["genera"]
    conda: 
        "envs/prokka.yaml"
    threads: 4
    log:
        "results/{genera}/logs/prokka/{sample}_prokka.log"
    shell:
        """
        prokka --force --cpus {threads} --kingdom Bacteria --genus {params.target} \
        --outdir "results/{params.genera}/prokka/{wildcards.sample}" --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.inputseqs} &> {log}
        """

rule roary:
    """
    Obtain core gene alignment and gene presence absence file
    """
    input:
        gff_files=expand("results/{genera}/prokka/{sample}/{sample}.gff", sample=SAMPLES, genera=config["genera"])
    output:
        "results/{genera}/roary/core_gene_alignment.aln",
        "results/{genera}/roary/gene_presence_absence.csv"
    conda:
        "envs/roary.yaml"
    params:
       outdir="results/{genera}/roary", 
       genera=config["genera"]
    threads: 4
    log:
        "results/{genera}/logs/roary/roary.log"
    shell:
        """
        rm -r results/{params.genera}/roary
        roary {input.gff_files} -e -n -v -p {threads} -f {params.outdir}
        """

rule fastree:
    """
    Generate newick file for Parnas
    """
    input:
        aln="results/{genera}/roary/core_gene_alignment.aln"
    output:
        nwk="results/{genera}/roary/tree.newick"
    conda:
        "envs/roary.yaml"
    log:
        "results/{genera}/logs/FastTree/FastTree.log"
    shell:
        """
        FastTree -nt -gtr {input.aln} > {output.nwk}
        """

rule fastani:
    """
    Run FastANI between our RefSeq and our queryseqs 
    """
    input:
        txt=config["txt"],
        ref = config["ref"]
    output:
        "results/{genera}/fastani/fastani.out"
    conda:
        "envs/fastani.yaml"
    log:
        "results/{genera}/logs/fastani/fastani.log"
    shell:
        """
        fastANI --ql {input.txt} -r {input.ref} -o {output} >> {log} 2>&1
        """