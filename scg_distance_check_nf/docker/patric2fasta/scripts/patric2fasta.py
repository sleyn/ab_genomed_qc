#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import argparse

# Create one sequence record
def make_record(record_id, record_seq):
    return SeqRecord(Seq(record_seq), id=record_id, description="")

# Main body
def main():
    parser = argparse.ArgumentParser(
            description="Convert PATRIC table to FASTA file")
    parser.add_argument('-t', '--table', type=str, action='store', help='PATRIC table')
    parser.add_argument('-o', '--out', type=str, action='store', help='Output FASTA file')
    parser.add_argument('--ref_sequence', type=str, action='store', default='-', help='If specified adds reference sequence to the FASTA file. Will be named as "reference"')
    args = parser.parse_args()
    
    input_table = args.table
    output_fasta = args.out
    
    patric_table = pd.read_csv(input_table, sep='\t', index_col=False)
    patric_table = patric_table[['feature.patric_id', 'feature.aa_sequence']]
    record_list = []
    
    patric_table.apply(
        lambda protein: record_list.append(
            make_record(protein['feature.patric_id'], protein['feature.aa_sequence'])
        ),
        axis=1
    )
    
    if args.ref_sequence != '-':
        record_list.append(
        	make_record("reference", args.ref_sequence)
        )
    
    with open(output_fasta, 'w') as output_file:
        SeqIO.write(record_list, output_file, "fasta")

if __name__ == "__main__":
    main()
