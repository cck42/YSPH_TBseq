### STATS ###



rule cleanstats:
    input: 
        flagstats = 'results/align/{sample}_{target}.flagstats',
        hoststats = 'results/clean/{sample}_{host}.flagstats',
        cleanstats = 'results/clean/{sample}_{host}_cleaned.flagstats',
    output:
        stats='results/stats/{sample}_{target}_{host}.alignstats',
    group: "statsgroup"
    params:
        mindepth=10,
        runtime=30,
    run:
        sample = wildcards.sample
        target = wildcards.target
        host = wildcards.host

        #get reads aligned
        reads = -1
        aligned = -1
        paired = -1
        with open(input.flagstats, "r") as f:
            for l in f:
                l = l.strip().split('\t')
                passreads = l[0]
                failreads = l[1]
                stat = l[2]
                #if stat == 'total (QC-passed reads + QC-failed reads)': 
                if stat == 'primary': 
                    reads = int(passreads)
                #elif stat == "mapped": 
                elif stat == "primary mapped": 
                    aligned = int(passreads)
                elif stat == "properly paired": 
                    paired = int(passreads)
        f.close()

        #get reads host aligned
        hreads = -1
        haligned = -1
        hpaired = -1
        with open(input.hoststats, "r") as f:
            for l in f:
                l = l.strip().split('\t')
                passreads = l[0]
                failreads = l[1]
                stat = l[2]
                #if stat == 'total (QC-passed reads + QC-failed reads)': 
                if stat == 'primary': 
                    hreads = int(passreads)
                #elif stat == "mapped": 
                elif stat == "primary mapped": 
                    haligned = int(passreads)
                elif stat == "properly paired": 
                    hpaired = int(passreads)
        f.close()

        #get cleaned reads
        creads = -1
        caligned = -1
        cpaired = -1
        with open(input.cleanstats, "r") as f:
            for l in f:
                l = l.strip().split('\t')
                passreads = l[0]
                failreads = l[1]
                stat = l[2]
                #if stat == 'total (QC-passed reads + QC-failed reads)': 
                if stat == 'primary': 
                    creads = int(passreads)
                #elif stat == "mapped": 
                elif stat == "primary mapped": 
                    caligned = int(passreads)
                elif stat == "properly paired": 
                    cpaired = int(passreads)
        f.close()

        #subprocess.run(["samtools","index",output.subsamp])
        f = open(output.stats, "w")
        print("\t".join(map(str,[sample, target,  reads,  aligned,   paired, round(aligned/reads,2), 
                                           host, hreads, haligned,  hpaired, round(hpaired/reads,2), 
                                        "clean", creads, caligned,  cpaired, round(creads/reads,2), 
                                        ])),file=f)
        f.close()


rule catstats:
    input:
        flagstats = expand(os.path.join(config['outdir'],'stats/{sample}_{{target}}_{{host}}.alignstats'),sample=SAMPLES)
    output:
        alstats=   'results/stats/{target}_{host}_alignstats.txt',
    group: "catstats"
    resources:
        mem_mb=8000,
        runtime=180,
    log:
        stderr="logs/stats/{target}_{host}_alignstats.err"
    shell:
        """
        cat {input.flagstats} > {output.alstats}
        """
