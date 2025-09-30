#/usr/bin/python3

import os
import sys
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Predict amplicons from primer pairs and a reference genome.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file with primer pairs.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome in FASTA format.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file for predicted amplicons.")
    return parser.parse_args()


