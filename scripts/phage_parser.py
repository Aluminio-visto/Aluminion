# %%
import pandas as pd
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# %%
def read_summary(infile):
       df = pd.read_table(infile, skipinitialspace=True, skiprows=32, sep=' ')
       df.drop(index=0, inplace=True)
       if len(df) > 0:
              df['MOST_COMMON_PHAGE_NAME(hit_genes_count)'] = df['MOST_COMMON_PHAGE_NAME(hit_genes_count)'].str.split(',').str[0].replace('\(.*\)', '', regex=True)
              df = df[['REGION', 'REGION_POSITION', 'REGION_LENGTH', 'COMPLETENESS(score)', 'SPECIFIC_KEYWORD',
                     'TOTAL_PROTEIN_NUM','PHAGE+HYPO_PROTEIN_PERCENTAGE',
                     'ATT_SITE_SHOWUP', 'MOST_COMMON_PHAGE_NAME(hit_genes_count)']]
              return df
       return None

# %%
def execute_blastn(assembly_fasta, phage_fna):
    mkbl_cmd = ['makeblastdb', '-in', assembly_fasta, '-parse_seqids', '-dbtype', 'nucl']
    subprocess.run(mkbl_cmd)
    blast_cmd = ['blastn', '-query', phage_fna, '-db', assembly_fasta, '-outfmt', '6 qseqid sseqid pident qcovhsp length qlen slen qstart qend sstart send sframe evalue bitscore']
    pipe = subprocess.Popen(blast_cmd, stdout=subprocess.PIPE)
    df_blast = pd.read_table(pipe.stdout, header=None)
    df_blast.columns = ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'sframe', 'evalue', 'bitscore']
    return df_blast

# %%
def process_blastn(df_blast):
    fagos_id = df_blast['qseqid'].unique()
    subdf = pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'sstart', 'send'])
    for fago_id in fagos_id:
        best_id = df_blast.loc[df_blast['qseqid'] == fago_id][['bitscore']].idxmax()
        selected_row = df_blast.loc[best_id, ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'sstart', 'send']]
        subdf = pd.concat([subdf, selected_row])
    subdf['qseqid'] =  subdf['qseqid'].astype(str)
    return subdf

def extract_fasta(df, phage_fna, output_dir):
    fasta_sequences = {record.id: record for record in SeqIO.parse(phage_fna, "fasta")}

    # Recorrer cada fila del DataFrame
    for _, row in df.iterrows():
        fago = row['Fago']
        cluster = row['Cluster']
        
        # Obtener la secuencia correspondiente al contig
        if fago not in fasta_sequences:
            print(f"Fago {fago} no encontrado en {phage_fna}.")
            continue
        
        full_sequence = fasta_sequences[fago].seq
        
        # Crear un registro SeqRecord
        record_id = f"{cluster}_{fago}_{sample}"
        seq_record = SeqRecord(
            full_sequence,
            id=record_id,
            description=f"{cluster}_{fago}_{sample}"
        )
        
        # Guardar el archivo FASTA
        output_path = f"{output_dir}/{record_id}.fasta"
        SeqIO.write(seq_record, output_path, "fasta")
        print(f"Archivo guardado: {output_path}")

    return None

# %%
if __name__ == "__main__":
    original_path = os.path.abspath(sys.argv[1])
    summary_df = pd.DataFrame(columns=['sample', 'Fago', 'contig', 'Start', 'End', 'length', 'Cluster', 'COMPLETENESS(score)',
                                        'SPECIFIC_KEYWORD', 'TOTAL_PROTEIN_NUM',
                                        'PHAGE+HYPO_PROTEIN_PERCENTAGE', 'ATT_SITE_SHOWUP'])
    for phage_path in glob.glob(f'{original_path}/09_phages/phastest_deep/*/'):
        # Extraer archivos necesarios
        sample = phage_path.split('/')[-2]
        phage_sum = os.path.join(phage_path, 'summary.txt')
        phage_fna = os.path.join(phage_path, 'region_DNA.txt')
        assembly_fasta = os.path.join(original_path, f'03_assemblies/{sample}/assembly.fasta')
        print(sample)
        # si tenemos todos los archivos disponibles
        if os.path.isfile(phage_sum) and os.path.isfile(phage_fna) and os.path.isfile(assembly_fasta):
            # procesar el summary
            df_phastest = read_summary(phage_sum)
            # si hay resultados, extraer la posición real de los fagos con un blastn frente al assembly original y sacar fastas de los fagos
            if df_phastest.empty == False:
                df_blast = execute_blastn(assembly_fasta, phage_fna)
                df_blast2 = process_blastn(df_blast)
                combo_df = df_phastest.merge(df_blast2, left_on='REGION', right_on='qseqid', how='outer')
                combo_df['sample'] = sample
                combo_df = combo_df[['sample', 'REGION', 'sseqid', 'sstart', 'send', 'length', 'MOST_COMMON_PHAGE_NAME(hit_genes_count)',
                                    'COMPLETENESS(score)', 'SPECIFIC_KEYWORD', 'TOTAL_PROTEIN_NUM', 'PHAGE+HYPO_PROTEIN_PERCENTAGE', 'ATT_SITE_SHOWUP']]
                combo_df.rename(columns={'REGION': 'Fago', 'sseqid': 'contig', 'sstart': 'Start', 'send': 'End', 'MOST_COMMON_PHAGE_NAME(hit_genes_count)': 'Cluster'}, inplace=True)
                extract_fasta(combo_df, phage_fna, f'{original_path}/09_phages/')
                summary_df = pd.concat([summary_df, combo_df], axis=0)
        else:
            continue

    summary_df.to_csv(f'{original_path}/09_phages/phage_summary.csv', index=False)
    print(f'Phage summary in {original_path}/09_phages/phage_summary.csv.')
