#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import subprocess
import shutil
from Bio import SeqIO
from BCBio import GFF
import os
import glob
import sys
import re
import warnings
from Bio import BiopythonDeprecationWarning
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

def abr_parse(abr_out):
    abr_raw = pd.read_table(abr_out)
    if len(abr_raw) > 0:
        abr_raw[['pos_beg', 'pos_end']] = pd.DataFrame(abr_raw['SEQUENCE'].str.split('_').str[-2:].tolist(), index=abr_raw.index)
        abr_raw.drop(abr_raw.loc[abr_raw['%IDENTITY'] < 90].index, inplace=True)
        abr_raw.drop(abr_raw.loc[abr_raw['%COVERAGE'] < 80].index, inplace=True)
        df_abr = abr_raw[['pos_beg', 'pos_end', 'GENE']].copy()  # .copy() to avoid SettingWithCopyWarning
        df_abr.columns = ['pos_beg', 'pos_end', 'abr_ann']
        df_abr['pos_beg'] = df_abr['pos_beg'].astype('int64')
        df_abr['pos_end'] = df_abr['pos_end'].astype('int64')
        return df_abr
    return pd.DataFrame(columns=['pos_beg', 'pos_end', 'abr_ann'])

def prokka_parse(prokka_dir):
    # Accumulate rows in a plain list and build the DataFrame once at the end.
    # The previous `df.loc[-1] = ...; df.index += 1` pattern silently overwrites
    # a row when the index value -1 already exists in the frame, which can happen
    # if pandas internals decide to renumber after concat.
    rows = []
    for prokka_file in glob.glob(f'{prokka_dir}/*.gff'):
        with open(prokka_file) as in_handle:
            for rec in GFF.parse(in_handle):
                start, end = rec.id.split('_')[-2:]
                for feature in rec.features:
                    gene = feature.qualifiers.get('gene', ['NA'])[0]
                    product = feature.qualifiers.get('product', ['NA'])[0]
                    rows.append([start, end, f'{gene};{product}'])
    df_prokka = pd.DataFrame(rows, columns=['pos_beg', 'pos_end', 'prokka_ann'])
    df_prokka['pos_beg'] = df_prokka['pos_beg'].astype('int64')
    df_prokka['pos_end'] = df_prokka['pos_end'].astype('int64')
    return df_prokka

def write_fasta(cds_output_file, d_int):
    with open(cds_output_file, "w") as cds_output_handle:
        for k, v in d_int.items():
            cds_output_handle.write(k + "\n")
            cds_output_handle.write(v + "\n")
    return None

def annotate_cds(cds_output_file, input_path, replicon, integron):
    prokka_dir = input_path + f'/prokka_{replicon}_{integron}'
    abr_path = input_path + f'/abricate_{replicon}_{integron}.out'
    # check=True so a Prokka or Abricate failure halts integron annotation for this
    # sample with a visible traceback instead of producing empty downstream files.
    prokka_cmd = ['prokka', '--quiet', '--force', cds_output_file, '--outdir', prokka_dir]
    subprocess.run(prokka_cmd, check=True)
    with open(abr_path, 'w') as abr_out:
        subprocess.run(['abricate', cds_output_file], stdout=abr_out, check=True)
    return prokka_dir, abr_path

def extract_fastas(gbk_file, cds_output_file, fna_file, integron):
    # Open GBK file and output files
    with open(gbk_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            flag = 0
            d_int = {}
            d_nucl = {}
            for feature in record.features:
                # Get start and end coordinates
                start = feature.location.start + 1  # +1 for 1-based coordinates
                end = feature.location.end
                # Extract the integron sequence
                if "integron" in feature.type and feature.qualifiers.get('integron_id')[0] == integron:
                    flag = 1
                    integron_seq = feature.extract(record.seq)
                    integron_header = f'>{integron}_{start}_{end}'
                    d_nucl[integron_header] = str(integron_seq)
                # Extract and save CDS sequences
                elif feature.type == "CDS" and flag == 1:
                    # Extract CDS sequence
                    cds_seq = feature.extract(record.seq)
                    # Get protein_id or locus_tag
                    protein_id = feature.qualifiers.get('protein_id', ['Unknown'])[0]
                    # Write to FASTA file
                    header = f">{protein_id}_{start}_{end}"
                    d_int[header] = str(cds_seq)

                elif "integron" in feature.type:
                    flag = 0

        # Integron nucleotide sequence (fna)
        write_fasta(fna_file, d_nucl)
        # Integron CDS sequences (faa)
        write_fasta(cds_output_file, d_int)
            
    return d_int

def extract_info(sample, subdf, replicon, integron, input_path, original_path):
    # Integrase
    try:
        integrase_row = subdf[subdf['annotation'] == 'intI'].iloc[0]
        integrase_model = integrase_row['model']
        integrase_strand = integrase_row['strand']
    except:
        integrase_model = ''
        integrase_strand = -1

    # attC sites
    try:
        attc_models = subdf[subdf['type_elt'] == 'attC']['model'].tolist()
    except:
        attc_models = ''

    try:
        # Extract faa and fna for integrons
        gbk_file = input_path + f'/{replicon}.gbk'
        cds_output_file = input_path + f'/{replicon}_{integron}.faa'
        fna_file = input_path + f'/{replicon}_{integron}.fna'
        d_int = extract_fastas(gbk_file, cds_output_file, fna_file, integron)
    except:
        print(f'\033[93m[WARNING]\033[0m GBK parsing error. Check the integrity of {gbk_file}')
        return None

    df_abr = pd.DataFrame(columns = ['pos_beg', 'pos_end', 'abr_ann'])
    df_prokka = pd.DataFrame(columns = ['pos_beg', 'pos_end', 'prokka_ann'])
    # If the integron has CDS, annotate them with Prokka and Abricate
    if len(d_int.keys()) > 0:
        prokka_dir, abr_out = annotate_cds(cds_output_file, input_path, replicon, integron)
        df_prokka = prokka_parse(prokka_dir)
        df_abr = abr_parse(abr_out)
    # Merge with annotations, prioritise Abricate and reorient df if needed
    subdf['pos_beg'] = subdf['pos_beg'].astype('int64')
    subdf['pos_end'] = subdf['pos_end'].astype('int64')
    mid_df = pd.merge(subdf, df_prokka, on=['pos_beg', 'pos_end'], how='outer')
    final_df = pd.merge(mid_df, df_abr, on=['pos_beg', 'pos_end'], how='outer')
    final_df['ann'] = final_df['abr_ann'].fillna(final_df['prokka_ann'])
    final_df.drop_duplicates(inplace=True)
    if integrase_strand > 0:
        final_df = final_df[::-1]
    
    # Group cassettes
    cassettes = []
    current_cassette = []

    for _, row in final_df.iterrows():
        if row['type_elt'] == 'attC':
            cassettes.append(current_cassette)
            current_cassette = []
        elif row['type_elt'] == 'protein' and row['annotation'] != 'intI':
            if pd.isnull(row['ann']):
                ann = "NA"  # Drop the "hypothetical protein" description
            else:
                ann = row['ann'].split(';')[0]  # Keep only the gene name
            current_cassette.append(ann)

    # In case the last attC is not recognised
    cassettes.append(current_cassette)
    # Remove empty cassettes
    cassettes = [i for i in cassettes if i != []]

    # Names already cleaned, build the gene string
    genes = '_'.join([i[0] for i in cassettes])
    genes = re.sub(r'[^a-zA-Z0-9\-\_]', '', genes)

    # Convert Python lists to clean comma-separated strings
    formatted_cassettes = [", ".join(c) for c in cassettes]

    # Data list for final output row
    contig = final_df['ID_replicon'].values[0]
    name = f'{integron}_{replicon}_{genes}_{sample}'
    start = min(final_df['pos_beg'])
    end = max(final_df['pos_end'])
    size = end - start
    type = final_df['type'].values[0]
    info = [sample, contig, name, size, start, end, type, integrase_model]

    # Append formatted cassettes
    for i in formatted_cassettes: info.append(i)
    info.extend([""] * (20-len(info)))

    # Save nucleotide sequence. shutil.copy is portable across OSes and raises
    # immediately on permission / missing-source errors instead of swallowing them.
    shutil.copy(fna_file, f'{original_path}/11_integrons/{name}.fasta')

    return info

def run_parsing(original_path, out_folder=None):
    """Main function called by the orchestrator script (parser.py)."""
    if out_folder is None:
        out_folder = original_path

    original_path = os.path.abspath(original_path)
    summary_df = pd.DataFrame(columns=['Sample', 'Pl/Chr', 'Name', 'Size', 'Start', 'End', 'Type', 'Integrase', 'Cassette 1',
                                        'Cassette 2', 'Cassette 3', 'Cassette 4', 'Cassette 5', 'Cassette 6',
                                        'Cassette 7', 'Cassette 8', 'Cassette 9', 'Cassette 10', 'Cassette 11', 'Cassette 12'])
    
    search_pattern = os.path.join(original_path, '11_integrons', '*', 'Results_Integron_Finder_*', '*.integrons')
    integron_files = glob.glob(search_pattern)

    if not integron_files:
        print(f"\033[93m[WARNING]\033[0m No integron files found to parse in: 11_integrons/")
        return

    for integron_file in integron_files:
        input_path = os.path.dirname(os.path.abspath(integron_file))
        # The sample directory is two levels up from the .integrons file:
        # 11_integrons/<sample>/Results_Integron_Finder_*/<file>.integrons
        # os.path keeps this portable across Windows and POSIX path separators.
        sample = os.path.basename(os.path.dirname(input_path))
        
        try:
            df_integron = pd.read_table(integron_file, comment='#')
        except:
            print(f'No integrons in {sample}')
            continue

        # Divide into integrons and chromosomes
        grouped = df_integron.groupby(['ID_replicon', 'ID_integron'])
        subdfs = {}

        for (replicon, integron), group in grouped:
            key = f"{replicon}_{integron}"
            subdfs[key] = group
            info = extract_info(sample, subdfs[key], replicon, integron, input_path, original_path)
            if info:
                summary_df.loc[len(summary_df)] = info

    # Save to the main folder (where the reporter will look for it)
    output_csv = os.path.join(out_folder, 'integron_summary.csv')
    summary_df.to_csv(output_csv, index=False)
    print(f" -> Integrons parsed successfully to {output_csv}")


if __name__ == "__main__":
    # Allows running this script independently from the command line
    if len(sys.argv) > 1:
        run_parsing(sys.argv[1])
    else:
        print("Usage: python integron_parser.py <run_directory>")