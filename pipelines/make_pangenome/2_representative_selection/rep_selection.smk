configfile: "config/SP_parnas.yaml"
RADIUSES = set(config['radius'])

rule all:
    input:
        expand("results/{genera}/parnas/reduced_tree.newick", genera=config["genera"]),
        expand("results/{genera}/parnas/representatives.txt", genera=config["genera"]), 
        expand("results/{genera}/parnas/{rad}/colors_tree.csv", rad = RADIUSES, genera=config["genera"]),
        expand("results/{genera}/parnas/{rad}/{rad}_cluster_reps.tab", rad = RADIUSES, genera=config["genera"]), 
        expand("results/{genera}/parnas/{rad}/colors_formatted.csv", rad = RADIUSES, genera=config["genera"]),
        expand("results/{genera}/parnas/{rad}/{rad}_subtree.tre", rad = RADIUSES, genera=config["genera"])

rule reduce_tree:
    """
    Remove refseq and outgroup from tree for rep selection
    """
    input:
        tree=config["tree"]
    output:
        reduced_tree="results/{genera}/parnas/reduced_tree.newick"
    params:
        ref=config["ref"],
        outgroup=config["outgroup"]
    shell:
        r"""
        Rscript scripts/reduce_tree.R {input.tree} {output.reduced_tree} {params.ref} {params.outgroup}
        """

rule estimate_representative_number:
    """
    This estimates how many representatives is sufficient to maximally cover diversity
    """
    input:
        rules.reduce_tree.output.reduced_tree
    output:
        scores = "results/{genera}/parnas/estimated_number_representatives.csv",
        representatives = "results/{genera}/parnas/representatives.txt",
        colortree = "results/{genera}/parnas/representative_strains.tre", 
        subtree1 = "results/{genera}/parnas/subtree.tre" 
    params: 
        max_num = config["max_num"]
    shell:
        """
        parnas -t {input} -n {params.max_num} \
        --diversity {output.scores} > {output.representatives} \
        --color {output.colortree} \
        --subtree {output.subtree1} 
        """

rule create_diversity_figure:
    """
    This creates a figure using the diversity covered data in the prior step
    """
    input:
        scores = rules.estimate_representative_number.output.scores
    output:
        figure = "results/{genera}/parnas/Estimated_Reps_Needed.jpg"
    shell:
        r"""
        Rscript "scripts/estimate_reps_needed.R" --input {input.scores} --output {output.figure}
        """

rule colortrees_with_radius:
    """
    Find a minimum number of representatives that cover most of the diversity across our tree
    """
    input:
        rules.reduce_tree.output.reduced_tree
    params:
        radius = "{rad}"
    output:
        colors = "results/{genera}/parnas/{rad}/colors_tree.csv",
        reps = "results/{genera}/parnas/{rad}/{rad}_cluster_reps.tab", 
        subtree2 = "results/{genera}/parnas/{rad}/{rad}_subtree.tre"
    shell:
        """
        parnas -t {input} --cover --radius {params.radius} \
        --subtree {output.subtree2} \
        --clusters {output.colors} > {output.reps} \
        """

rule combine_colortrees:
    """
    Combine reps with numerous radiuses from prior rule 
    """
    input:
        colors = rules.colortrees_with_radius.output.colors, 
        reps = rules.colortrees_with_radius.output.reps
    output:
        formatted = "results/{genera}/parnas/{rad}/colors_formatted.csv"
    shell:
        r"""
        Rscript scripts/colorcombiner.R \
        --representatives-file {input.reps} \
        --colors-file {input.colors} \
        --output {output}
        """