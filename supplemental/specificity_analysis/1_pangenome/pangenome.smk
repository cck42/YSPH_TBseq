import pandas as pd

configfile: "config/TB_Specificity.yaml"

samples_df = pd.read_csv("tsv/TB_Specificity.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
SEQ = {row.sample_id: {"sample_id": row.sample_id, "seq": row.seq_path} for row in samples_df.itertuples()} 

rule all:
    input:
        expand("results/{genera}/specificity_analysis/pangenome/genome_annotation/{sample}/{sample}.gff", genera=config["genera"],sample=SAMPLES),
        expand("results/{genera}/specificity_analysis/pangenome/pangenome_calc/core_gene_alignment.aln", genera=config["genera"]),
        expand("results/{genera}/specificity_analysis/pangenome/pangenome_calc/gene_presence_absence.csv", genera=config["genera"]),
        expand("results/{genera}/specificity_analysis/pangenome/pangenome_calc/tree.newick", genera=config["genera"]), 
        expand("results/{genera}/specificity_analysis/pangenome/average_nucleotide_identity/fastani.out", genera=config["genera"])

rule genome_annotation:
    """
    Annotate genomes and obtain gff files required to run Roary
    """
    input:
        inputseqs=lambda wildcards: SEQ[wildcards.sample]["seq"]
    output:
        "results/{genera}/specificity_analysis/pangenome/genome_annotation/{sample}/{sample}.gff"
    params:
        target=config["target"], 
        genera=config["genera"],
        kingdom="Bacteria",
        outdir="results/{genera}/specificity_analysis/pangenome/genome_annotation/{sample}"
    threads: 4
    log:
        stdout = "logs/{genera}/specificity_analysis/pangenome/genome_annotation/{sample}/prokkka.out",
        stderr = "logs/{genera}/specificity_analysis/pangenome/genome_annotation/{sample}/prokka.err"
    shell:
        """
        module unload miniconda 
        source activate /home/flg9/.conda/envs/prokka

        prokka \
        --force --cpus {threads} --kingdom {params.kingdom} --genus {params.target} \
        --outdir {params.outdir} --prefix {wildcards.sample} \
        --locustag {wildcards.sample} {input.inputseqs} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule pangenome_calc:
    """
    Calculate pangenome content to obtain core gene alignment 
    and gene presence absence file that wher'll use for plotting later
    """
    input:
        gff_files=expand("results/{genera}/specificity_analysis/pangenome/genome_annotation/{sample}/{sample}.gff", sample=SAMPLES, genera=config["genera"])
    output:
        "results/{genera}/specificity_analysis/pangenome/pangenome_calc/core_gene_alignment.aln",
        "results/{genera}/specificity_analysis/pangenome/pangenome_calc/gene_presence_absence.csv"
    params:
       outdir="results/{genera}/specificity_analysis/pangenome/pangenome_calc", 
       genera=config["genera"]
    threads: 4
    log:
        stdout = "logs/{genera}/specificity_analysis/pangenome/pangenome_calc/roary.out",
        stderr = "logs/{genera}/specificity_analysis/pangenome/pangenome_calc/roary.err"
    shell:
        """
        module unload miniconda 
        source activate /home/flg9/.conda/envs/roary

        # Roary has a bug where it creates the output dir and then errors
        # out because the outdir already exists so we are bypassing that here

        rm -r {params.outdir}

        roary \
        {input.gff_files} -e -n -v -p {threads} -f {params.outdir} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule ml_tree:
    """
    Generate newick file for Parnas
    """
    input:
        aln="results/{genera}/specificity_analysis/pangenome/pangenome_calc/core_gene_alignment.aln"
    output:
        nwk="results/{genera}/specificity_analysis/pangenome/pangenome_calc/tree.newick"
    log:
        stdout = "logs/{genera}/specificity_analysis/pangenome/pangenome_calc/fasttree.out",
        stderr = "logs/{genera}/specificity_analysis/pangenome/pangenome_calc/fasttree.err"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/roary

        FastTree \
        -nt -gtr {input.aln} > {output.nwk} \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule average_nucleotide_identity:
    """
    Run FastANI between our RefSeq and our queryseqs to get ANI
    """
    input:
        txt=config["txt"],
        ref = config["ref"]
    output:
        "results/{genera}/specificity_analysis/pangenome/average_nucleotide_identity/fastani.out"
    log:
        stdout = "logs/{genera}/specificity_analysis/pangenome/average_nucleotide_identity/fastani.out",
        stderr = "logs/{genera}/specificity_analysis/pangenome/average_nucleotide_identity/fastani.err"
    shell:
        """
        module unload miniconda 
        source activate /vast/palmer/pi/turner/flg9/conda_envs/fastani

        fastANI \
        --ql {input.txt} -r {input.ref} -o {output} \
        1>> {log.stdout} 2>> {log.stderr}
        """