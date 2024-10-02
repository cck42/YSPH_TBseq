

ruleorder: mykrobe_skeleton > mykrobe

rule mykrobe_skeleton:
    input:
        R1 = expand('results/rawdata/{sample}_R1.fastq.gz',sample=SAMPLE1),
        R2 = expand('results/rawdata/{sample}_R2.fastq.gz',sample=SAMPLE1),
    params:
        outprefix=expand(os.path.join(config['outdir'],'mykrobe/{sample}_mykrobe'),sample=SAMPLE1),
        sample=SAMPLE1
    output:
        #mykrobepredict_json=os.path.join(config['outdir'],'mykrobe/{sample}_mykrobe.json'),
        #mykrobepredict_csv=os.path.join(config['outdir'],'mykrobe/{sample}_mykrobe.csv'),
        skeleton=directory(os.path.join(config['outdir'],'mykrobe_skeleton'))
    log:
        stdout="logs/mykrobe/skeleton.out",
        stderr="logs/mykrobe/skeleton.err"
    container: "docker://sethnr/pgcoe_bacseq:0.02"
    resources:
        disk_mb=10000,
        mem_mb=16000,
        runtime=120,
        cores=1
    shell:
        """
        mykrobe predict \
        --sample {params.sample} \
        --species tb \
        --output {params.outprefix} \
        --format json_and_csv \
        --skeleton_dir ./results/mykrobe_skeleton \
        --tmp ./results/mykrobe/ \
        --seq {input.R1} {input.R2} 1> {log.stdout} 2> {log.stderr}
        """



rule mykrobe:
    input:
        R1 = 'results/rawdata/{sample}_R1.fastq.gz',
        R2 = 'results/rawdata/{sample}_R2.fastq.gz',
        skeleton=os.path.join(config['outdir'],'mykrobe_skeleton')
    params:
        outprefix=os.path.join(config['outdir'],'mykrobe/{sample}_mykrobe'),
    output:
        mykrobepredict_json=os.path.join(config['outdir'],'mykrobe/{sample}_mykrobe.json'),
        mykrobepredict_csv=os.path.join(config['outdir'],'mykrobe/{sample}_mykrobe.csv')
    log:
        stdout="logs/mykrobe/{sample}.out",
        stderr="logs/mykrobe/{sample}.err"
    container: "docker://sethnr/pgcoe_bacseq:0.02"
    resources:
        disk_mb=10000,
        mem_mb=16000,
        runtime=120,
        cores=1
    shell:
        """
        mykrobe predict \
        --sample {wildcards.sample} \
        --species tb \
        --output {params.outprefix} \
        --format json_and_csv \
        --skeleton_dir ./results/mykrobe_skeleton \
        --tmp ./results/mykrobe/ \
        --seq {input.R1} {input.R2} 1> {log.stdout} 2> {log.stderr}
        """


rule mykrobe_combine:
    input:
        mykrobepredict=expand(config['outdir']+'/mykrobe/{sample}_mykrobe.csv',sample=SAMPLES)
    params:
    output:
        mykrobe_combined=os.path.join(config['outdir'],'mykrobe/combined.csv'),
        fin="results/mykrobe/complete"
    log:
        stdout="logs/mykrobe/combine.out",
        stderr="logs/mykrobe/combine.err"
    shell:
        """
        echo -e '"sample","drug","susceptibility","variants (dna_variant-AA_variant:ref_kmer_count:alt_kmer_count:conf) [use --format json for more info]","genes (prot_mut-ref_mut:percent_covg:depth) [use --format json for more info]","mykrobe_version","files","probe_sets","genotype_model","kmer_size","phylo_group","species","lineage","phylo_group_per_covg","species_per_covg","lineage_per_covg","phylo_group_depth","species_depth","lineage_depth"' > {output.mykrobe_combined}
        for file in {input.mykrobepredict}; do
            tail -n +2 "$file" >> {output.mykrobe_combined}
        done 2> {log.stderr}


        touch {output.fin}
        """