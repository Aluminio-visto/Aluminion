#!/usr/bin/env python3

# Jorge R Grande - HUMV - Santander
# IS_parser.py takes an IS.tsv table from blastn
# and outputs a parsed table without repeated elements

import os
import pandas as pd
import numpy as np
import argparse


def get_arguments():
    parser = argparse.ArgumentParser(prog='IS_parser.py', description='IS_parser.py is part of the Aluminion analysis pipeline.')

    input_group = parser.add_argument_group('Input', 'Input parameters')
    input_group.add_argument('-i', '--input_folder', dest="input_folder", required=True,
                             help="Required. Input folder containing IS_chr.tsv", type=os.path.abspath)

    output_group = parser.add_argument_group('Output', 'Output parameters')
    output_group.add_argument('-o', '--out_dir', dest='out_dir', required=False,
                              help='Output folder (defaults to input_folder)', type=os.path.abspath)

    return parser.parse_args()


def main():
    args = get_arguments()

    input_folder = args.input_folder
    out_folder   = args.out_dir if args.out_dir else input_folder

    report     = os.path.join(input_folder, "IS_chr.tsv")
    report_out = os.path.join(out_folder,   "IS_chr_out.tsv")

    df = pd.read_csv(report, sep='\t', names=['IS', 'contig', 'start', 'end', '%ID', 'mismatch', 'evalue'])
    df['start2'] = df['start']
    df['end2']   = df['end']
    df = df.round({'start2': -2})
    df = df.round({'end2':   -2})
    df = df.sort_values(by=["start2", "mismatch"], ascending=True)
    df = df.drop_duplicates(subset='start2', keep='first')
    df = df.sort_values(by=["end2", "mismatch"], ascending=True)
    df = df.drop_duplicates(subset='end2', keep='first')
    df = df[df['%ID'] > 90]
    df = df.drop(columns=['contig', 'evalue', 'start2', 'end2'])

    df.to_csv(report_out, sep='\t')


if __name__ == '__main__':
    main()
