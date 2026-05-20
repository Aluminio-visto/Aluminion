#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def read_summary(infile):
    """Reads Phastest summary output.

    The header row is located by searching for the line that starts with 'REGION'
    instead of relying on a fixed skiprows offset — Phastest has shifted the
    number of preamble lines between releases, and a hardcoded value silently
    misaligns every column when the format changes.
    """
    header_idx = None
    with open(infile) as fh:
        for i, line in enumerate(fh):
            if line.lstrip().startswith('REGION '):
                header_idx = i
                break
    if header_idx is None:
        return None

    df = pd.read_table(infile, skipinitialspace=True, skiprows=header_idx, sep=r'\s+', engine='python')
    # Phastest writes a divider line of dashes right under the header; drop it.
    if len(df) > 0 and str(df.iloc[0, 0]).startswith('-'):
        df = df.drop(index=df.index[0]).reset_index(drop=True)

    if len(df) > 0:
        df['MOST_COMMON_PHAGE_NAME(hit_genes_count)'] = df['MOST_COMMON_PHAGE_NAME(hit_genes_count)'].str.split(',').str[0].replace(r'\(.*\)', '', regex=True)
        df = df[['REGION', 'REGION_POSITION', 'REGION_LENGTH', 'COMPLETENESS(score)', 'SPECIFIC_KEYWORD',
                'TOTAL_PROTEIN_NUM','PHAGE+HYPO_PROTEIN_PERCENTAGE',
                'ATT_SITE_SHOWUP', 'MOST_COMMON_PHAGE_NAME(hit_genes_count)']]
        return df
    return None

def execute_blastn(assembly_fasta, phage_fna):
    """Executes BLASTN to find exact phage regions in the original assembly."""
    mkbl_cmd = ['makeblastdb', '-in', assembly_fasta, '-parse_seqids', '-dbtype', 'nucl']
    # check=True surfaces makeblastdb failures (e.g. truncated FASTA, missing perms)
    # instead of letting the subsequent blastn fail with a less helpful error.
    subprocess.run(mkbl_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True)
    blast_cmd = ['blastn', '-query', phage_fna, '-db', assembly_fasta, '-outfmt', '6 qseqid sseqid pident qcovhsp length qlen slen qstart qend sstart send sframe evalue bitscore']
    pipe = subprocess.Popen(blast_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    df_blast = pd.read_table(pipe.stdout, header=None)
    pipe.wait()
    if pipe.returncode != 0:
        err = pipe.stderr.read().decode(errors='replace')
        raise RuntimeError(f"blastn failed (rc={pipe.returncode}): {err}")
    if not df_blast.empty:
        df_blast.columns = ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'sframe', 'evalue', 'bitscore']
    return df_blast

def process_blastn(df_blast):
    """Processes BLASTN output to extract the best hits."""
    if df_blast.empty:
        return pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'sstart', 'send'])
        
    fagos_id = df_blast['qseqid'].unique()
    subdf = pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'sstart', 'send'])
    for fago_id in fagos_id:
        best_id = df_blast.loc[df_blast['qseqid'] == fago_id]['bitscore'].idxmax()
        selected_row = df_blast.loc[[best_id], ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'sstart', 'send']]
        subdf = pd.concat([subdf, selected_row], ignore_index=True)
    subdf['qseqid'] = subdf['qseqid'].astype(str)
    return subdf

def extract_fasta(df, phage_fna, output_dir, sample):
    """Extracts the FASTA sequences of the identified phages."""
    fasta_sequences = {record.id: record for record in SeqIO.parse(phage_fna, "fasta")}

    for _, row in df.iterrows():
        fago = row['Fago']
        cluster = row['Cluster']
        
        if fago not in fasta_sequences:
            print(f" [WARNING] Phage {fago} not found in {phage_fna}.")
            continue
        
        full_sequence = fasta_sequences[fago].seq
        record_id = f"{cluster}_{fago}_{sample}"
        seq_record = SeqRecord(
            full_sequence,
            id=record_id,
            description=f"{cluster}_{fago}_{sample}"
        )
        
        output_path = os.path.join(output_dir, f"{record_id}.fasta")
        SeqIO.write(seq_record, output_path, "fasta")

    return None

def run_parsing(original_path, out_folder=None):
    """Main function called by the orchestrator (parser.py)."""
    if out_folder is None:
        out_folder = original_path

    original_path = os.path.abspath(original_path)
    summary_df = pd.DataFrame(columns=['Sample', 'Fago', 'contig', 'Start', 'End', 'length', 'Cluster', 'COMPLETENESS(score)',
                                        'SPECIFIC_KEYWORD', 'TOTAL_PROTEIN_NUM',
                                        'PHAGE+HYPO_PROTEIN_PERCENTAGE', 'ATT_SITE_SHOWUP'])
    
    search_pattern = os.path.join(original_path, '09_phages', 'phastest_deep', '*')
    
    for phage_path in glob.glob(search_pattern):
        if not os.path.isdir(phage_path):
            continue
            
        sample = os.path.basename(os.path.normpath(phage_path))
        phage_sum = os.path.join(phage_path, 'summary.txt')
        phage_fna = os.path.join(phage_path, 'region_DNA.txt')
        assembly_fasta = os.path.join(original_path, f'03_assemblies/{sample}/assembly.fasta')
        
        if os.path.isfile(phage_sum) and os.path.isfile(phage_fna) and os.path.isfile(assembly_fasta):
            df_phastest = read_summary(phage_sum)
            if df_phastest is not None and not df_phastest.empty:
                df_blast = execute_blastn(assembly_fasta, phage_fna)
                df_blast2 = process_blastn(df_blast)
                
                if not df_blast2.empty:
                    combo_df = df_phastest.merge(df_blast2, left_on='REGION', right_on='qseqid', how='outer')
                    combo_df['Sample'] = sample
                    combo_df = combo_df[['Sample', 'REGION', 'sseqid', 'sstart', 'send', 'length', 'MOST_COMMON_PHAGE_NAME(hit_genes_count)',
                                        'COMPLETENESS(score)', 'SPECIFIC_KEYWORD', 'TOTAL_PROTEIN_NUM', 'PHAGE+HYPO_PROTEIN_PERCENTAGE', 'ATT_SITE_SHOWUP']]
                    combo_df.rename(columns={'REGION': 'Fago', 'sseqid': 'contig', 'sstart': 'Start', 'send': 'End', 'MOST_COMMON_PHAGE_NAME(hit_genes_count)': 'Cluster'}, inplace=True)
                    
                    extract_fasta(combo_df, phage_fna, os.path.join(original_path, '09_phages'), sample)
                    summary_df = pd.concat([summary_df, combo_df], axis=0, ignore_index=True)

    output_csv = os.path.join(out_folder, 'phage_summary.csv')
    summary_df.to_csv(output_csv, index=False)
    print(f" -> Parsed phages saved to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        run_parsing(sys.argv[1])
    else:
        print("Usage: python phage_parser.py <run_directory>")