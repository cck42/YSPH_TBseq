import argparse
from Bio import SeqIO
import pandas as pd 
from Bio.Seq import Seq

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('pos_bed')
parser.add_argument('neg_bed')
parser.add_argument('bedfasta_pos')
parser.add_argument('bedfasta_neg')
args = parser.parse_args()

# load primer bed files 
pos_primers = pd.read_csv(args.pos_bed, sep='\t', header=None)
neg_primers = pd.read_csv(args.neg_bed, sep='\t', header=None)

# extract sequences from bed 
sequences_pos = []
sequences_neg = []
for df, sequences in [(pos_primers, sequences_pos), (neg_primers, sequences_neg)]:
    for _, row in df.iterrows():
        seq = Seq(row[6]) 
        name = '_'.join(row.astype(str))
        sequences.append(SeqIO.SeqRecord(seq, id=name, description=""))

# write sequences to a fasta file 
SeqIO.write(sequences_pos, args.bedfasta_pos, "fasta")
SeqIO.write(sequences_neg, args.bedfasta_neg, "fasta")