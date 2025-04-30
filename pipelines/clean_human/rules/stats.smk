### STATS ###



rule alignstats:
    input: 
        flagstats = 'results/align/{sample}_{target}.flagstats',
    output:
        stats='results/align/{sample}_{target}_alignstats.txt',
    group: "statsgroup"
    params:
        mindepth=10,
        runtime=30,
    run:
        sample = wildcards.sample
        target = wildcards.target

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
                if stat == 'total (QC-passed reads + QC-failed reads)': 
                    reads = int(passreads)
                elif stat == "mapped": 
                    aligned = int(passreads)
                elif stat == "properly paired": 
                    paired = int(passreads)
        f.close()


        #subprocess.run(["samtools","index",output.subsamp])
        f = open(output.stats, "w")
        print("\t".join(map(str,[sample, target, reads, aligned, paired])),file=f)
        f.close()


