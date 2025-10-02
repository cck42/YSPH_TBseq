#/usr/bin/python3

import os
import sys
import argparse
import pandas as pd
import pysam
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description="Predict amplicons from primer pairs and a reference genome.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file with primer pairs.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome in FASTA format.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file for predicted amplicons.")
    return parser.parse_args()


def parse_minimap_output(minimap_file):
    """
    Parses minimap2 output in PAF format and returns a DataFrame.
    """
    columns = [
        "query_name", "query_length", "query_start", "query_end",
        "strand", "target_name", "target_length", "target_start", "target_end",
        "residue_matches", "alignment_block_length", "mapping_quality"
    ]
    records = []
    with open(minimap_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            record = dict(zip(columns, fields[:12]))
            records.append(record)
    df = pd.DataFrame(records)
    df[['query_length', 'query_start', 'query_end', 'target_length', 'target_start', 'target_end',
        'residue_matches', 'alignment_block_length', 'mapping_quality']] = df[[
            'query_length', 'query_start', 'query_end', 'target_length', 'target_start', 'target_end',
            'residue_matches', 'alignment_block_length', 'mapping_quality']].apply(pd.to_numeric)
    return df  

def predict_amplicons(primer_fasta, reference_file, output_file, MAXLEN=3000):
    # Step 1: Run minimap2 to align primers to the reference genome
    minimap_output = f'{output_file}.paf'
    #minimap kmer 5, window 3, match 4, mismatch 2, min chain score 5
    minimap_command = [
        "minimap2", "-cxsr", "-k6", "-w4", "-A4", "-B2", "-m2",
        reference_file, primer_fasta
    ]
    print(f"Running command: {' '.join(minimap_command)} > {minimap_output}", file=sys.stderr)
    with open(minimap_output, "w") as out_f:
        subprocess.run(minimap_command, stdout=out_f, check=True)

    # Step 2: Parse minimap2 output
    alignments = parse_minimap_output(minimap_output)
    
    # Step 3: Identify potential amplicons
    amplicons = []
    fwd_aligns = alignments[alignments['strand'] == "+"]
    rev_aligns = alignments[alignments['strand'] == "-"]
        
    for _, fwd in fwd_aligns.iterrows():
        for _, rev in rev_aligns.iterrows():
            if fwd['target_end'] < rev['target_start']:
                amplicon_length = rev['target_start'] - fwd['target_end']
                if amplicon_length > MAXLEN:
                    continue
                else:
                    amplicons.append({
                        "p1": fwd['query_name'],
                        "p2": rev['query_name'],
                        "p1mismatches": fwd['query_length']-fwd['residue_matches'],
                        "p2mismatches": rev['query_length']-rev['residue_matches'],
                        "p1identity": round((fwd['residue_matches']/fwd['alignment_block_length']),2),
                        "p2identity": round((rev['residue_matches']/rev['alignment_block_length']),2),
                        #"p1alignlen": fwd['alignment_block_length'],
                        #"p2alignlen": rev['alignment_block_length'],
                        "target_name": fwd['target_name'],
                        "amplicon_start": fwd['target_end'],
                        "amplicon_end": rev['target_start'],
                        "amplicon_length": amplicon_length
                    })
    amplicon_df = pd.DataFrame(amplicons)
    amplicon_df.to_csv(output_file, sep="\t", index=False)

def primers_bed_to_fasta(bed_file, fasta_odd, fasta_even):
    """
    Reads primers from a BED file (sequence in column 5) and writes two FASTA files:
    one for odd-numbered primer pairs, one for even-numbered pairs.
    """
    odd = []
    even = []
    with open(bed_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            primer_name = fields[3]
            sequence = fields[6]
            idx = int(primer_name.split('_')[-2]) #name format = primername_idx_left/right

            fasta_entry = f">{primer_name}\n{sequence}\n"
            if (idx % 2) == 0:
                odd.append(fasta_entry)
            else:
                even.append(fasta_entry)
    with open(fasta_odd, "w") as f:
        f.writelines(odd)
    with open(fasta_even, "w") as f:
        f.writelines(even)


def estimate_genome_coverage(amplicon1, amplicon2, reference_fasta, maxmismatch=3):
    """
    Estimates total genome coverage from predicted amplicons.
    Returns the fraction of the reference genome covered by at least one amplicon.
    """

    # Load reference genome lengths
    ref_lengths = {}
    with pysam.FastaFile(reference_fasta) as ref:
        for name in ref.references:
            ref_lengths[name] = ref.get_reference_length(name)
    total_length = sum(ref_lengths.values())

    # Load amplicons
    amplicons1 = pd.read_csv(amplicon1, sep="\t")
    amplicons2 = pd.read_csv(amplicon2, sep="\t")
    covered = {name: set() for name in ref_lengths}

    for _, row in amplicons1.iterrows():
        chrom = row['target_name']
        start = int(row['amplicon_start'])
        end = int(row['amplicon_end'])
        if(row['p1mismatches'] > maxmismatch or row['p2mismatches'] > maxmismatch):
            continue
        if chrom in covered:
            covered[chrom].update(range(start, end))

    for _, row in amplicons2.iterrows():
        chrom = row['target_name']
        start = int(row['amplicon_start'])
        end = int(row['amplicon_end'])
        if(row['p1mismatches'] > maxmismatch or row['p2mismatches'] > maxmismatch):
            continue
        if chrom in covered:
            covered[chrom].update(range(start, end))

    total_covered = sum(len(bases) for bases in covered.values())
    coverage_fraction = total_covered / total_length if total_length > 0 else 0
    return (total_covered, coverage_fraction)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict amplicons from primer pairs and a reference genome.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file with primer pairs.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome in FASTA format.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file for predicted amplicons.")
    args = parser.parse_args()

    primers_bed_to_fasta(args.input, "p1.fasta", "p2.fasta")

    predict_amplicons("p1.fasta", args.reference, f'{args.output}_p1', MAXLEN=4000)
    predict_amplicons("p2.fasta", args.reference, f'{args.output}_p2', MAXLEN=4000)

    total_covered, coverage = estimate_genome_coverage(f'{args.output}_p1', f'{args.output}_p2', args.reference)

    with open(f"{args.output}_coverage.txt", "w") as cov_out:
        print(f"Bases covered: {total_covered}, Estimated {args.reference} coverage by amplicons: {coverage:.4f}", file=sys.stderr)
        print(f"Bases covered: {total_covered}, Estimated {args.reference} coverage by amplicons: {coverage:.4f}", file=cov_out)

    # Clean up temporary files
    os.remove("minimap_output.paf")
    os.remove("p1.fasta")
    os.remove("p2.fasta")