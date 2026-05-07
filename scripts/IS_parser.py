#!/home/usuario/miniconda3/envs/nano/bin/python

#############################################################
# Jorge R Grande - HUMV - Santander
# IS_parser.py takes an IS.tsv table from blastn
# and outputs a parsed table without repeated elements



import os
import pandas as pd
import numpy as np
import argparse


def get_arguments():

    """
    Esta función permite introducir los parámetros al llamar al programa en la terminal
    """

    parser = argparse.ArgumentParser(prog = 'parser.py', description = 'Parser.py is part of the ION analysis pipeline.')

    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input_folder', dest="input_folder", required=True, help="Required. Input folder to be summarized", type=os.path.abspath)

    # input_group.add_argument('-db', '--pubMLST_database', dest="pubMLST_database", required=False, default="/home/usuario/miniconda3/envs/mlst/db/pubmlst/", help="Folder containing the MLST schemas", type=os.path.abspath)

    output_group = parser.add_argument_group('Output', 'Output parameters')

    output_group.add_argument('-o', '--out_dir', dest='out_dir', required=False, help='Final output folder for reports, same as input_folder by default', type=os.path.abspath)

    arguments = parser.parse_args()

    return arguments

args = get_arguments()


##################
# Taking in args:#
##################

input_folder        = args.input_folder

if args.out_dir:
    out_folder      = args.out_folder
else:
    out_folder      = input_folder


report              = input_folder + "/IS_chr.tsv"

report_out          = out_folder + "/IS_chr_out.tsv"

##################
# Taking in args:#
##################


df = pd.read_csv(report, sep='\t', names=['IS','contig','start','end','%ID','mismatch','evalue'])
df['start2'] = df['start']
df['end2']  = df['end']
df = df.round({'start2':-2})
df = df.round({'end2' :-2})
df = df.sort_values(by=["start2","mismatch"], ascending=True)
df = df.drop_duplicates(subset='start2', keep='first')
df = df.sort_values(by=["end2","mismatch"], ascending=True)
df = df.drop_duplicates(subset='end2', keep='first')
df = df[df['%ID'] > 90]
df = df.drop(columns=['contig','evalue','start2','end2'])

df.to_csv(report_out, sep='\t')