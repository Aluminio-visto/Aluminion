# %%
import pandas as pd
from glob import glob
import dateutil.parser
import re
import ast
import os
import argparse
import numpy as np

# %%
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to process sequencing runs and check reads and assembly QC data")
    
    # Required argument for 'base_run'
    parser.add_argument("--input_path", type=str, help="Path to the sequencing run directory")

    # Optional argument for 'output_run'
    parser.add_argument("--output_file", type=str, default="data_seq_new.tsv",
                        help="Output filename (default: data_seq_new.tsv)")
    parser.add_argument("--init", action='store_true',
                        help="First run: create data_seq.tsv and data_analysis.tsv from scratch "
                             "(also triggered automatically when those files do not exist)")
    parser.add_argument("--alert-all-mge", action='store_true',
                        help="Alert on any shared MGE regardless of clinical relevance. "
                             "Default: only alert on plasmids or integrons carrying resistance genes.")

    return parser.parse_args()

# %%
def parse_minion_sum(summary):
    d_sum = {}
    with open(summary, "r") as f:
        for line in f:
            if line.startswith("instrument"):
                d_sum['instrument'] = line.split(sep='=')[-1].rstrip('\n')
            if line.startswith("flow_cell_id"):
                d_sum['flow_cell'] = line.split(sep='=')[-1].rstrip('\n')
            if line.startswith("protocol="):
                d_sum['flow_cell_type'] = line.split(sep=':')[0].split(sep='_', maxsplit=1)[1].rstrip('\n')
                d_sum['barcoding_kit'] = line.split(sep=':')[-2].rstrip('\n')
            if line.startswith("started"):
                d_sum['start'] = dateutil.parser.parse(line.split(sep='=')[-1])
            if line.startswith("acquisition_stopped"):
                d_sum['end'] = dateutil.parser.parse(line.split(sep='=')[-1])

        duration = d_sum['end'] - d_sum['start']                       # For build-in functions
        d_sum['fecha'] = d_sum['start'].strftime('%Y-%m-%d')  # sequencing date
        hours, minutes, seconds = str(duration).split(sep=':')
        d_sum['duracion'] = hours + 'h ' + minutes + 'min.'  # run duration

    return d_sum

# %%
def parse_minion_report(report):
    d_report = {}
    with open(report, "r") as rep:
        doc= rep.read()
        items = re.findall(r'\"total_pores\":\"\d*\"', doc)
        d_report['poros_ini'] = items[0].split(":")[-1].strip("\"")
        d_report['poros_fin'] = items[-1].split(":")[-1].strip("\"")
        
    return d_report

# %%
def get_gfa(gfa):
    with open(gfa, 'r') as file:

        d_links = {}
        d_paths = {}

        for line in file.readlines():
            # Parse links
            if line.startswith('L'):
                edge1 = line.split()[1]
                edge2 = line.split()[3]
                if edge1 not in d_links:
                    d_links[edge1] = []
                d_links[edge1].append(edge2)
                d_links[edge1] = list(set(d_links[edge1]))
            # Parse paths
            if line.startswith('P'):
                contig = line.split()[1]
                edges = line.split()[2].split(',')
                edges = [re.sub(r"[+-]", "", edge) for edge in edges]
                d_paths[contig] = set(edges)
    return d_links, d_paths

# %%
def parse_info(file_path, d_links, d_paths):
    df = pd.read_table(file_path)
    df.sort_values('length', ascending=False, inplace=True)    # Default order, but just in case
    d_info = {}

    # Is the largest contig closed?
    d_info['contig1_closed'] = df['circ.'][0] == 'Y'

    # Size ratio (2 largest contigs / total)
    contig1 = df['#seq_name'][0]
    if len(df) != 1:
        total_length = sum(df['length'])
        l_contig1 = df['length'][0]
        l_contig2 = df['length'][1]
        d_info['ratio_total'] = (l_contig1 + l_contig2) / total_length

        # Ratio between the two largest contigs
        d_info['ratio_greatest_contigs'] = l_contig2 / l_contig1

        # Check if the two largest contigs share a common node
        contig2 = df['#seq_name'][1]
        edges1 = d_paths[contig1]
        edges2 = d_paths[contig2]
        d_info['are_linked'] = False
        for k1 in edges1:
            for k2 in d_links.get(k1, ['no_edge']):
                if any(k2 in e2 for e2 in edges2):
                    d_info['are_linked'] = True
    else:
        d_info['ratio_total'] = 1
        d_info['ratio_greatest_contigs'] = 0
        d_info['are_linked'] = False

    # Number of plasmids
    if d_info['are_linked']:
        chromosomes = {contig1, contig2}
    else:
        chromosomes = {contig1}
    plasmids = [cont for cont in df['#seq_name'] if cont not in chromosomes]
    d_info['n_plas'] = len(plasmids)

    # Number of closed plasmids
    closed_plas = df[df['#seq_name'].isin(plasmids)]['circ.']
    closed_plas.reset_index()
    d_info['n_closed_plas'] = sum(closed_plas == 'Y')

    # return n_contigs, n_closed, contig1_closed, ratio_total, ratio_greatest_contigs, are_linked, contig1, contig2, n_plas, closed_plas
    return d_info

# %%
# Assembly score calculation
def calculate_score(d_info):

    score = 0
    d_score = {}

    # Score for chromosome quality
    # Ratio of the two largest contigs to total length
    if d_info['ratio_total'] > 0.8:
        score += 1.5
        d_score['ratio'] = 1.5
    else:
        d_score['ratio'] = 0

    # If the second longest contig is at most 10% of the largest
    if d_info['ratio_greatest_contigs'] < 0.1:
        score += 1.25
        d_score['ratio_1_2'] = 1.25
        # If the chromosome is closed
        if d_info['contig1_closed']:
            score += 0.75
            d_score['contig1_closed'] = 0.75
        else:
            d_score['contig1_closed'] = 0
    # Otherwise, check if contigs are linked (insertion sequence, 0.5) or not
    else:
        d_score['ratio_1_2'] = 0
        if d_info['are_linked'] == True:
            score += 0.75
            d_score['greatest_linked'] = 0.75
        else:
            d_score['greatest_linked'] = 0

    # Score for plasmids
    plasmid_score = 1.5
    if d_info['n_plas'] > 0:
        plasmid_score = d_info['n_closed_plas'] * 1.5 / d_info['n_plas']
    score += plasmid_score
    d_score['plasmid_score'] = plasmid_score

    score = round(score, 2)

    return score, d_score

# %%
def get_assembly_score(lista_cepas, base_run):
    d_quality = {}
    for sample in lista_cepas['ID']:
        assembly_path = f'{base_run}/03_assemblies/{sample}/'
        info = os.path.join(assembly_path, 'assembly_info.txt')
        gfa = os.path.join(assembly_path, 'assembly_graph.gfa')

        # If assembly exists, calculate score; otherwise return 0
        if os.path.isfile(info):
            d_links, d_paths = get_gfa(gfa)
            d_info = parse_info(info, d_links, d_paths)
            score, d_score = calculate_score(d_info)
        else:
            score = 0
            d_score = {}
        d_quality[sample] = score
    
    return d_quality

# %%
def _is_amr_gene(name):
    """Return True if a normalised cassette gene name looks like an AMR gene.

    Excludes: emrE (efflux pump without clinical AMR significance), hypothetical
    proteins logged as 'na', and IS-element transposases.
    Everything else (bla*, aac*, aph*, dfr*, sul*, tet*, cat*, qnr*, aad*, …)
    is treated as an AMR gene.
    """
    n = name.lower()
    return n not in ('emre', 'na', '') and not n.startswith('is')


def _parse_cassette_cell(cell):
    """Extract base gene names from an integron cassette cell (Python list string)."""
    if not cell or (isinstance(cell, float) and pd.isna(cell)):
        return []
    try:
        genes = ast.literal_eval(str(cell))
        if not isinstance(genes, list):
            genes = [str(genes)]
    except (ValueError, SyntaxError):
        genes = [str(cell)]
    result = []
    for gene in genes:
        name = gene.split(';')[0].strip()
        name = re.sub(r'_\d+$', '', name)
        if name and name.upper() != 'NA':
            result.append(name.lower())
    return result


def build_mge_table(df_copla, df_integron, df_phage, seq_date, alert_all=False):
    """Build a unified long-format MGE table for one sequencing run.

    Produces one row per MGE occurrence per sample. Used to detect
    MGEs shared across sequential runs.

    alert_all: if False (default), only plasmids/integrons carrying resistance
    genes are included. If True, any mobilisable plasmid or cassette-carrying
    integron is included regardless of resistance content.
    """
    rows = []
    trivial = {'-', '', 'nan'}

    # Plasmids: one row per contig, using a composite identity key.
    # Key hierarchy: PTU (copla classification) > Rep+MOB+MPF combination.
    # Plasmids where MOB, MPF and AbR are all trivial are always skipped —
    # no meaningful key can be formed and there is nothing clinically relevant.
    # Without alert_all, plasmids without resistance genes are also skipped.
    if df_copla is not None and not df_copla.empty:
        for _, row in df_copla.iterrows():
            mob = str(row.get('MOB', '-')).strip()
            mpf = str(row.get('MPF', '-')).strip()
            abr = str(row.get('AbR', '-')).strip()

            if mob in trivial and mpf in trivial and abr in trivial:
                continue  # no key possible and no clinical signal — always skip
            if not alert_all and abr in trivial:
                continue  # no resistance genes — skip unless alert_all

            ptu = str(row.get('PTU', '-')).strip()
            rep = str(row.get('Rep', '-')).strip()

            if ptu and ptu not in trivial:
                mge_key = ptu
            else:
                parts = []
                if rep not in trivial:
                    rep_sorted = '|'.join(
                        sorted(r.strip() for r in rep.split(';')
                               if r.strip() not in trivial))
                    if rep_sorted:
                        parts.append(f'Rep:{rep_sorted}')
                if mob not in trivial:
                    parts.append(f'MOB:{mob}')
                if mpf not in trivial:
                    parts.append(f'MPF:{mpf}')
                mge_key = ';'.join(parts)
                if not mge_key:
                    continue

            rows.append({
                'ID': row['Sample'],
                'Seq_date': seq_date,
                'MGE_type': 'plasmid',
                'MGE_key': mge_key,
                'MGE_detail': abr if abr not in trivial else '',
            })

    # Integrons: one row per complete integron with at least one cassette.
    # Without alert_all, integrons with no AMR genes in their cassettes are skipped.
    cassette_cols = [f'Cassette {i}' for i in range(1, 13)]
    if df_integron is not None and not df_integron.empty:
        for _, row in df_integron.iterrows():
            if str(row.get('Type', '')) != 'complete':
                continue
            genes = []
            for col in cassette_cols:
                if col in row:
                    genes.extend(_parse_cassette_cell(row[col]))
            if not genes:
                continue
            if not alert_all and not any(_is_amr_gene(g) for g in genes):
                continue  # only housekeeping genes (emrE, etc.) — skip
            rows.append({
                'ID': row['Sample'],
                'Seq_date': seq_date,
                'MGE_type': 'integron',
                'MGE_key': '|'.join(sorted(genes)),
                'MGE_detail': '',
            })

    # Prophages: one row per intact phage (matched by cluster accession)
    if df_phage is not None and not df_phage.empty:
        for _, row in df_phage.iterrows():
            completeness = str(row.get('COMPLETENESS(score)', ''))
            cluster = str(row.get('Cluster', ''))
            if completeness.startswith('intact') and cluster and cluster != 'nan':
                rows.append({
                    'ID': row['Sample'],
                    'Seq_date': seq_date,
                    'MGE_type': 'prophage',
                    'MGE_key': cluster,
                    'MGE_detail': '',
                })

    cols = ['ID', 'Seq_date', 'MGE_type', 'MGE_key', 'MGE_detail']
    return pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)


def find_shared_mges(current_mge, history_mge):
    """Return MGEs in the current run that also appear in historical data.

    Matching is exact on (MGE_type, MGE_key). Self-matches (same sample ID
    in both current and history, e.g. a repeated sequencing) are excluded.
    """
    empty = pd.DataFrame(columns=[
        'ID_current', 'ID_historical', 'Seq_date_historical',
        'MGE_type', 'MGE_key', 'MGE_detail',
    ])
    if history_mge.empty or current_mge.empty:
        return empty

    hist = history_mge[['ID', 'Seq_date', 'MGE_type', 'MGE_key']].rename(
        columns={'ID': 'ID_historical', 'Seq_date': 'Seq_date_historical'})
    curr = current_mge[['ID', 'MGE_type', 'MGE_key', 'MGE_detail']].rename(
        columns={'ID': 'ID_current'})

    merged = pd.merge(curr, hist, on=['MGE_type', 'MGE_key'])
    shared = merged[merged['ID_current'] != merged['ID_historical']].drop_duplicates()
    return shared.sort_values(['MGE_type', 'MGE_key']).reset_index(drop=True)


# %%

def main():
    # Input path con todos los archivos
    args = parse_arguments()
    base_run = args.input_path
    output_run = args.output_file

    # Column schemas for the two output databases
    COLS_SEQ = [
        "Lab_id", "Strain", "ID", "Barcode", "Barcode_rep1", "Barcode_rep2",
        "Seq_date", "Seq_date_rep1", "Seq_date_rep2", "DNA_conc",
        "Depth", "Assembly_score", "Extraction_kit", "Barcoding_kit",
        "Instrument", "Flowcell_type", "Flowcell", "Pores_start", "Pores_end",
        "Seq_hours", "Samples_per_run", "Samples_to_repeat",
        "Median_length_pre", "Median_quality_pre", "N_reads_pre", "N_bases_pre",
        "Median_length_post", "Median_quality_post", "N_reads_post", "N_bases_post",
        "Pct_bases_kept",
    ]
    COLS_ANALYSIS = [
        "Lab_id", "ID", "Barcode", "Depth", "Assembly_score",
        "Subspecies", "MLST", "Serotype", "KO_locus", "Contaminants",
        "Carbapenemase", "ESBL", "Other_resistance", "N_AMR_genes", "AMRscore", "VIRscore",
        "Plasmids", "Prophages", "Integrons",
        "allele_1", "allele_2", "allele_3", "allele_4", "allele_5", "allele_6", "allele_7",
        "Possible_MLSTs", "Possible_alleles", "Majority_genus", "Majority_species", "MLST_scheme",
    ]

    # %%
    # Input files
    # MinION final run summaries
    summary_matches = glob(os.path.join(base_run, "final_summary*txt"))
    report_matches  = glob(os.path.join(base_run, "report_*.json"))
    summary = summary_matches[0] if summary_matches else ""
    report  = report_matches[0]  if report_matches  else ""

    # Quality control statistics
    qc_r = os.path.join(base_run, "QC_reads.csv")
    qc_a = os.path.join(base_run, "QC_assembly.csv")

    # Historical sequencing database
    tabla = os.path.join(base_run, "data_seq.tsv")

    # Historical analysis database
    anali = os.path.join(base_run, "data_analysis.tsv")

    # Per-run sample list (distinct from data_seq.tsv — see list_seq.tsv template)
    cepas = os.path.join(base_run, "list_seq.tsv")

    # Taxonomy
    taxon = os.path.join(base_run, "taxonomy.csv")

    # EGMs
    plasmids   = os.path.join(base_run, "copla_modif.csv")
    fagos      = os.path.join(base_run, "phage_summary.csv")
    integrones = os.path.join(base_run, "integron_summary.csv")

    # Auto-detect first run: if historical files are absent, behave as --init
    init_mode = args.init or not os.path.isfile(tabla) or not os.path.isfile(anali)
    if init_mode:
        print("[INFO] First-run mode: historical database files not found or --init passed. "
              "Creating data_seq.tsv and data_analysis.tsv from scratch.")

    # Check that all required files exist
    required_files = [f for f in [summary, report, qc_r, qc_a, cepas, taxon] if f]
    missing_files = [f for f in required_files if not os.path.isfile(f)]
    if missing_files:
        raise FileNotFoundError(f"The following files do not exist: {missing_files}")


    # %%
    # Load inputs
    # Table with sequenced and pending samples
    if init_mode:
        datos_seq = pd.DataFrame(columns=COLS_SEQ)
    else:
        datos_seq = pd.read_csv(tabla, sep='\t')
        # Drop empty columns (Unnamed: )
        datos_seq = datos_seq.loc[:, ~datos_seq.columns.str.contains('^Unnamed: ')]
        # Set column data types
        datos_seq['Lab_id'] = datos_seq['Lab_id'].astype(str)
        datos_seq["Barcode"]    = datos_seq["Barcode"].astype(str)
        datos_seq['Barcode'] = datos_seq['Barcode'].replace('nan', np.nan)

    # Table with current run sample information
    lista_cepas = pd.read_csv(cepas, sep='\t', usecols=["Lab_id", "Strain", "ID", "Barcode", "DNA_conc"], dtype={'Barcode': 'string'})
    lista_cepas['Barcode'] = lista_cepas['Barcode'].str.replace(r'barcode', '', regex=True)

    # Add new samples
    nuevas_filas = lista_cepas[~lista_cepas['ID'].isin(datos_seq['ID'])]
    datos_seq = pd.concat([datos_seq, nuevas_filas], ignore_index=True)

    # Dicts with technical run metadata (optional: only present for complete runs)
    d_sum    = parse_minion_sum(summary)    if (summary    and os.path.isfile(summary))    else {}
    d_report = parse_minion_report(report)  if (report     and os.path.isfile(report))     else {}

    # Tables with sequencing and assembly quality
    QC_reads = pd.read_csv(qc_r, sep='\t', decimal='.', thousands=',')
    QC_assembly = pd.read_csv(qc_a, sep='\t')

    # Dictionary with assembly information
    d_quality = get_assembly_score(lista_cepas, base_run)

    # %%
    # Output files
    # Sequencing data output
    output_run = os.path.join(base_run, output_run)
    analisis_run = os.path.join(base_run, "data_analysis_new.tsv")

    # Initialize the output table with the current run samples and technical metadata
    columnas = ["Lab_id", "Strain", "ID", "Seq_date", "Seq_date_rep1", "Seq_date_rep2",
                "Extraction_kit", "Barcoding_kit", "Barcode", "Barcode_rep1", "Barcode_rep2", "Instrument",
                "Flowcell_type", "Flowcell", "Pores_start", "Pores_end", "Seq_hours", "Samples_per_run",
                "Samples_to_repeat", "Yield_Mbp", "is_repeated", "Temp_C", "Voltage", "Reads_per_hour",
                "Mbp_per_hour", "N50_kbp"]

    # Create 'result' by copying 'lista_cepas' and adding missing columns as null
    result = lista_cepas.copy()
    for col in columnas:
        if col not in result.columns:
            result[col] = pd.NA

    # %%
    # Populate table with technical metadata (NaN when MinION files are absent)
    result["Seq_date"]       = d_sum.get('fecha',          pd.NA)
    result["Barcoding_kit"]  = d_sum.get('barcoding_kit',  pd.NA)
    result["Extraction_kit"] = "DNeasy Blood & Tissue"
    result["Instrument"]     = d_sum.get('instrument',     pd.NA)
    result["Flowcell_type"]  = d_sum.get('flow_cell_type', pd.NA)
    result["Flowcell"]       = d_sum.get('flow_cell',      pd.NA)
    result["Seq_hours"]      = d_sum.get('duracion',       pd.NA)
    result["Pores_start"]    = d_report.get('poros_ini',   pd.NA)
    result["Pores_end"]      = d_report.get('poros_fin',   pd.NA)

    # %%
    # Populate with QC_reads.csv data
    QC_reads = QC_reads.rename(columns={"Sample":"ID",
                            "Median length" : "Median_length_pre",
                            "Median quality" : "Median_quality_pre",
                            "Total reads" : "N_reads_pre",
                            "Total bases" : "N_bases_pre",
                            "Median length.1" : "Median_length_post",
                            "Median quality.1" :"Median_quality_post",
                            "Total reads.1" : "N_reads_post",
                            "Total bases.1" : "N_bases_post"})

    QC_reads = QC_reads.drop(columns=["MaxQ", "Longest read", "Sample.1", "Samp", "MaxQ.1", "Longest read.1"], errors='ignore')
    QC_reads[["N_bases_post"]].apply(pd.to_numeric)
    result2 = pd.merge(result, QC_reads, on="ID", how='outer')

    # Populate with QC_assembly.csv data
    QC_assembly = QC_assembly.rename(columns={"Sample":"ID"})
    QC_assembly["ratio"] = QC_assembly["Largest contig"]/QC_assembly["Total length"]
    QC_assembly = QC_assembly.drop(columns=["GC (%)", "# predicted genes (>= 300 bp)"], errors='ignore')

    result3 = pd.merge(result2, QC_assembly, on="ID", how='outer')

    result3[["Total length", "N_bases_post"]].apply(pd.to_numeric)

    result3["Depth"] = result3["N_bases_post"].div(result3["Total length"])
    result3["Depth"] = result3["Depth"].round(0).astype('Int64')

    result3["Pct_bases_kept"] = result3["N_bases_post"].div(result3["N_bases_pre"])

    # Add assembly quality score
    result3['Assembly_score'] = result3['ID'].map(d_quality)



    orden_final = ["Lab_id", "Strain", "ID", "Barcode", "Barcode_rep1", "Barcode_rep2", "Seq_date", "Seq_date_rep1", "Seq_date_rep2", "DNA_conc",
            "Depth", 'Assembly_score', "Extraction_kit", "Barcoding_kit", "Instrument", "Flowcell_type",
            "Flowcell", "Pores_start", "Pores_end", "Seq_hours", "Samples_per_run", "Samples_to_repeat",
            "Median_length_pre", "Median_quality_pre", "N_reads_pre", "N_bases_pre",
            "Median_length_post", "Median_quality_post", "N_reads_post", "N_bases_post", "Pct_bases_kept"]
    result3 = result3[orden_final]

    # %%
    Ncepas_inicial = lista_cepas.shape[0]
    Ncepas_bien = (result3["Depth"] > 30.0).sum()
    Ncepas_repetir = Ncepas_inicial - Ncepas_bien

    result3["Samples_per_run"]      = Ncepas_inicial
    result3["Samples_to_repeat"] = Ncepas_repetir

    result3['Lab_id'] = result3['Lab_id'].astype(str)
    result3["Barcode"]    = result3["Barcode"].astype(str)

    # Merge tables on 'ID'
    merged_df = pd.merge(datos_seq, result3, on='ID', how='left', suffixes=('', '_result3'))

    # Assign values for "Seq_date_rep2" and "Barcode_rep2"
    merged_df['Seq_date_rep2'] = merged_df['Seq_date_rep2'].combine_first(
        merged_df.apply(lambda x: x['Seq_date_result3'] if pd.notna(x['Seq_date_rep1']) else None, axis=1))
    merged_df['Barcode_rep2'] = merged_df['Barcode_rep2'].combine_first(
        merged_df.apply(lambda x: x['Barcode_result3'] if pd.notna(x['Barcode_rep1']) else None, axis=1))

    # Assign values for "Seq_date" and "Barcode"
    merged_df['Seq_date_rep1'] = merged_df['Seq_date_rep1'].combine_first(
        merged_df.apply(lambda x: x['Seq_date_result3'] if pd.notna(x['Seq_date']) else None, axis=1))
    merged_df['Barcode_rep1'] = merged_df['Barcode_rep1'].combine_first(
        merged_df.apply(lambda x: x['Barcode_result3'] if pd.notna(x['Barcode']) else None, axis=1))

    # Fill any remaining missing values in the original columns
    merged_df['Seq_date'] = merged_df['Seq_date'].combine_first(merged_df['Seq_date_result3'])
    merged_df['Barcode'] = merged_df['Barcode'].combine_first(merged_df['Barcode_result3'])






    # %%
    # Fill empty rows in datos_seq with values from result3
    for column in datos_seq.columns:
        if column not in ['ID', 'Barcode', 'Barcode_rep1', 'Barcode_rep2', 'Seq_date', 'Seq_date_rep1', 'Seq_date_rep2']:  # Skip join key columns
            merged_df[column] = merged_df[column + '_result3'].combine_first(merged_df[column])

    # Elimina las columnas extra de result3
    merged_df = merged_df[datos_seq.columns]

    # Round integer columns to remove the ugly ".0"
    int_columns = [
        'Depth','Samples_per_run','Samples_to_repeat','Median_length_pre',
        'N_reads_pre','N_bases_pre','Median_length_post','N_reads_post','N_bases_post'
    ]
    for col in int_columns:
        if col in merged_df.columns:
            merged_df[col] = merged_df[col].astype(str).str.replace(',', '', regex=False).str.strip()
            merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce').round().astype('Int64')
    # Round the last column to 2 decimal places:
    merged_df["Pct_bases_kept"] = (merged_df["Pct_bases_kept"]*100).round(1)
    seq_out_path = os.path.join(base_run, "data_seq.tsv") if init_mode else output_run
    merged_df.to_csv(seq_out_path, index=False, sep="\t")
    print(f" -> data_seq written to: {seq_out_path}")

    # %%
    if not init_mode:
        analisis = pd.read_csv(anali, sep='\t')
        analisis.rename(columns={"Muestra":"Lab_id", "Serotipo": "Serotype"}, inplace=True)
    # (in init mode analisis is built from result4 alone — see below)

    taxon2 = pd.read_csv(taxon, sep= ',')
    taxon2.rename(columns={"Sample":"ID"}, inplace=True)
    taxon2.drop_duplicates(subset=['ID'], keep='first', inplace=True)

    result4 = pd.merge(taxon2, lista_cepas, on="ID", how='outer')
    result4.drop(columns=["Strain", "DNA_conc"], inplace=True)

    # Load MGE files (None when absent)
    df_pl    = pd.read_csv(plasmids,   sep=',') if os.path.isfile(plasmids)   else None
    df_fagos = pd.read_csv(fagos,      sep=',') if os.path.isfile(fagos)      else None
    df_int   = pd.read_csv(integrones, sep=',') if os.path.isfile(integrones) else None

    if df_fagos is not None:
        df_fagos.rename(columns={'sample': 'Sample'}, inplace=True)

    # Mobile Genetic Elements (MGEs) — count per sample for data_analysis
    result4[["Plasmids", "ICEs", "Prophages", "Integrons"]] = 0
    result4["ICEs"] = '0'  # not yet implemented

    if df_pl is not None:
        result4['Plasmids'] = result4['ID'].map(df_pl['Sample'].value_counts(), na_action='ignore')
    if df_fagos is not None:
        result4['Prophages'] = result4['ID'].map(df_fagos['Sample'].value_counts(), na_action='ignore')
    if df_int is not None:
        result4['Integrons'] = result4['ID'].map(df_int['Sample'].value_counts(), na_action='ignore')

    result4[['Plasmids', 'ICEs', 'Prophages', 'Integrons']] = result4[['Plasmids', 'ICEs', 'Prophages', 'Integrons']].fillna(0)

    result4 = result4.merge(merged_df[['ID', 'Assembly_score', 'Depth']], on='ID', how='inner')

    nwo = COLS_ANALYSIS
    result4 = result4.reindex(columns=nwo)
    result4['ID'] = result4['ID'].astype(str)
    result4["Barcode"]  = result4["Barcode"].astype(str)

    if init_mode:
        # First run: no historical data to merge with, result4 IS the database
        analisis_final = result4
    else:
        analisis = analisis[nwo]
        analisis['ID'] = analisis['ID'].astype(str)
        analisis["Barcode"]  = analisis["Barcode"].astype(str)
        analisis['Barcode']  = analisis['Barcode'].replace('nan', np.nan)

        # Merge historical data with new data
        analisis_final = pd.merge(analisis, result4, on='ID', how='left', suffixes=('', '_result4'))
        for column in analisis.columns:
            if column != 'ID':
                analisis_final[column] = analisis_final[column + '_result4'].combine_first(analisis_final[column])
        analisis_final = analisis_final[analisis.columns]

    # Replace "barcode13" with just "13"
    analisis_final['Barcode'] = analisis_final['Barcode'].str.replace(r'barcode', '', regex=True)

    # Round integer columns to remove the ugly ".0"
    int_columns = [
        "Lab_id", 'Depth', "N_AMR_genes", "AMRscore", "VIRscore",
        "Plasmids", "Prophages", "Integrons",
    ]
    for col in int_columns:
        if col in analisis_final.columns:
            analisis_final[col] = analisis_final[col].astype(str).str.replace(',', '', regex=False).str.strip()
            analisis_final[col] = pd.to_numeric(analisis_final[col], errors='coerce').round().astype('Int64')

    analysis_out_path = os.path.join(base_run, "data_analysis.tsv") if init_mode else analisis_run
    analisis_final.to_csv(analysis_out_path, index=False, sep='\t')
    print(f" -> data_analysis written to: {analysis_out_path}")

    # ── MGE cross-run comparison ──────────────────────────────────────────────
    mge_db_path = os.path.join(base_run, "data_mge.tsv")
    current_mge = build_mge_table(df_pl, df_int, df_fagos, d_sum.get('fecha', ''),
                                  alert_all=args.alert_all_mge)

    if init_mode:
        current_mge.to_csv(mge_db_path, index=False, sep='\t')
        print(f" -> MGE database created: {mge_db_path}")
    else:
        if os.path.isfile(mge_db_path):
            history_mge = pd.read_csv(mge_db_path, sep='\t')
            shared = find_shared_mges(current_mge, history_mge)
            if not shared.empty:
                shared_path = os.path.join(base_run, "mge_shared.tsv")
                shared.to_csv(shared_path, index=False, sep='\t')
                print(f" -> Shared MGEs: {len(shared)} events found → {shared_path}")
            else:
                print(" -> No shared MGEs found with historical runs")
            updated_mge = pd.concat([history_mge, current_mge], ignore_index=True)
        else:
            print("[WARN] data_mge.tsv not found; creating from this run only (no comparison possible)")
            updated_mge = current_mge
        updated_mge.to_csv(mge_db_path, index=False, sep='\t')
        print(f" -> MGE database updated: {mge_db_path}")

# %%
if __name__ == "__main__":
    main()
