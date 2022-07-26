# PATRIC to FASTA

Script to convert PATRIC table containing `feature.patric_id` and `feature.aa_sequence` columns to FASTA format.

Script can add a reference sequence with ID `reference`.

```
usage: patric2fasta.py [-h] [-t TABLE] [-o OUT] [--ref_sequence REF_SEQUENCE]

Convert PATRIC table to FASTA file

options:
  -h, --help            show this help message and exit
  -t TABLE, --table TABLE
                        PATRIC table
  -o OUT, --out OUT     Output FASTA file
  --ref_sequence REF_SEQUENCE
                        If specified adds reference sequence to the FASTA file. Will be named as "reference"
```