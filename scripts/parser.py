#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Aluminion — parser.py
# Author : Jorge R. Grande — HUMV / IDIVAL, Santander
# Purpose: Aggregate the per-tool outputs of aluminion.sh into a single set of
#          consolidated tables (taxonomy, AMR, MLST) plus an Excel summary.
# =============================================================================

import os
import argparse
import warnings

import pandas as pd

from _log import get_logger

log = get_logger(__name__)

# -----------------------------------------------------------------------------
# Optional sub-parser modules (imported lazily so a missing script doesn't crash
# the whole consolidation step — the affected sub-stage is just skipped).
# -----------------------------------------------------------------------------
try:
    import phage_parser
    import integron_parser
    import copla_parser
except ImportError as e:
    warnings.warn(f"Missing script in scripts/ folder: {e}")


# =============================================================================
# HELPERS
# =============================================================================
def safe_read_csv(filepath, required_cols, sep=',', **kwargs):
    """Read a CSV defensively.

    Returns an empty DataFrame with the requested ``required_cols`` if the file
    is missing, empty, or unparsable. Missing columns are added (as empty
    object columns) to the returned frame.
    """
    if not os.path.exists(filepath):
        log.warning('File not found, skipping: %s', filepath)
        return pd.DataFrame(columns=required_cols)
    try:
        if os.stat(filepath).st_size == 0:
            log.warning('Empty file: %s', filepath)
            return pd.DataFrame(columns=required_cols)

        df = pd.read_csv(filepath, sep=sep, **kwargs)
        for col in required_cols:
            if col not in df.columns:
                df[col] = pd.Series(dtype='object')
        return df
    except Exception as e:
        log.warning('Error reading %s: %s', filepath, e)
        return pd.DataFrame(columns=required_cols)


def get_arguments():
    parser = argparse.ArgumentParser(
        prog='parser.py',
        description='Aluminion Parser: aggregate per-tool outputs into a WGS summary.',
    )

    input_group = parser.add_argument_group('Input', 'Input parameters')
    input_group.add_argument('-i', '--input_folder', dest='input_folder', required=True,
                             help='Required. Folder containing aluminion.sh outputs to summarise.',
                             type=os.path.abspath)
    input_group.add_argument('-db', '--pubMLST_database', dest='pubMLST_database', required=False,
                             default='/home/usuario/miniconda3/envs/mlst/db/pubmlst/',
                             help='Folder containing PubMLST schemas.',
                             type=os.path.abspath)

    output_group = parser.add_argument_group('Output', 'Output parameters')
    output_group.add_argument('-o', '--out_dir', dest='out_dir', required=False,
                              help='Final output folder.', type=os.path.abspath)

    skip_group = parser.add_argument_group('Skip Modules', 'Bypass specific analysis parsers')
    skip_group.add_argument('--skip-phages', action='store_true',
                            help='Skip phage parsing (Phastest).')
    skip_group.add_argument('--skip-integrons', action='store_true',
                            help='Skip integron parsing (Integron Finder).')
    skip_group.add_argument('--skip-plasmids', action='store_true',
                            help='Skip plasmid parsing (Copla).')
    skip_group.add_argument('--skip-typing', action='store_true',
                            help='Skip typing data (MLST, Kleborate, ECTyper, GAMBIT).')
    skip_group.add_argument('--skip-kraken', action='store_true',
                            help='Skip Kraken2 contamination parsing.')
    skip_group.add_argument('--skip-abr', action='store_true',
                            help='Skip resistance gene extraction (Abricate).')

    return parser.parse_args()


def preflight_check(input_folder, args):
    """Warn the user about missing input files before any processing starts."""
    checks = []

    if not args.skip_abr:
        checks.append(('AbR_report.csv',
                       os.path.join(input_folder, 'AbR_report.csv'),
                       '--skip-abr'))
    if not args.skip_kraken:
        checks.append(('kraken2/species.csv',
                       os.path.join(input_folder, '04_taxonomies/kraken2/species.csv'),
                       '--skip-kraken'))
        checks.append(('kraken2/genus.csv',
                       os.path.join(input_folder, '04_taxonomies/kraken2/genus.csv'),
                       '--skip-kraken'))
    if not args.skip_typing:
        checks.append(('gambit.csv',
                       os.path.join(input_folder, '04_taxonomies/gambit.csv'),
                       '--skip-typing'))
        checks.append(('kleborate output',
                       os.path.join(input_folder,
                                    '04_taxonomies/kleborate/enterobacterales__species_output.txt'),
                       '--skip-typing'))
        checks.append(('ectyper output.tsv',
                       os.path.join(input_folder, '04_taxonomies/ectyper/output.tsv'),
                       '--skip-typing'))
        checks.append(('mlst.csv',
                       os.path.join(input_folder, 'mlst.csv'),
                       '--skip-typing'))

    missing = [(name, flag) for name, path, flag in checks if not os.path.exists(path)]
    if missing:
        log.warning('Preflight: the following input files were not found:')
        for name, flag in missing:
            log.warning('  ✗  %s  (skip with %s)', name, flag)
        log.warning('Preflight: processing will continue — missing modules will '
                    'produce empty tables.')
    else:
        log.info('Preflight: all expected input files found.')


# =============================================================================
# MLST PROCESSING
# =============================================================================
MLST_HEADER = [
    'Sample', 'MLST_scheme', 'MLST',
    'allele_1', 'allele_2', 'allele_3', 'allele_4',
    'allele_5', 'allele_6', 'allele_7', 'allele_8',
    'Possible_MLSTs', 'Possible_alleles',
]


def _process_mlst_row(row_fields, pubmlst_db_root):
    """Process a single line from ``mlst.csv``.

    Returns a (possibly extended) list. Padding to the full schema is the
    caller's responsibility, since rows with no MLST scheme legitimately have
    fewer than ``len(MLST_HEADER)`` fields.

    Three cases:

    1. No MLST scheme for the species (``MLST_scheme == '-'``): rewrite column 1
       and return the row as-is.
    2. Scheme present and ST resolved: return the row as-is.
    3. Scheme present but ST unresolved (``MLST == '-'``): look up the closest
       PubMLST allele combinations and append them in the last two columns.
    """
    # The first column from `mlst` is a file path — strip directories and the
    # trailing extension so the value matches sample IDs elsewhere.
    raw_path = row_fields[0]
    sample_id = raw_path.split('/')[1].split('.')[0] if '/' in raw_path else raw_path
    row_fields[0] = sample_id
    organism = row_fields[1]

    if organism == '-':
        row_fields[1] = 'No associated MLST scheme'
        return row_fields

    if row_fields[2] != '-':
        return row_fields

    # Scheme known but ST unresolved — lookup against the PubMLST schema.
    gene_to_allele = {}
    for cell in row_fields[3:]:
        if '(' not in cell:
            continue
        gene = cell.split('(')[0]
        allele = cell.split('(')[1].strip(')')
        if allele.isdigit():
            gene_to_allele[gene] = allele

    safe_genes = list(gene_to_allele.keys())
    db_path = os.path.join(pubmlst_db_root, organism, f'{organism}.txt')

    if safe_genes and os.path.exists(db_path):
        query = pd.DataFrame([gene_to_allele]).astype('int64')
        schema = pd.read_csv(db_path, sep='\t')
        candidates = pd.merge(query, schema, how='left', on=safe_genes)
    else:
        candidates = pd.DataFrame([{'ST': []}])

    # When the merge produced a single all-NaN row, take only the ST column;
    # otherwise take the column right after ST (the allele combination column).
    if candidates.shape[0] == 1 and candidates['ST'].isna().any():
        candidate_alleles = candidates.iloc[:, candidates.columns.get_indexer(['ST'])]
    else:
        candidate_alleles = candidates.iloc[:, candidates.columns.get_indexer(['ST']) + 1]

    possible_alleles = candidate_alleles.to_dict('list')
    for gene, alleles in possible_alleles.items():
        possible_alleles[gene] = [
            'No allele for this gene matches this combination' if pd.isna(a) else a
            for a in alleles
        ][:9]

    candidates['ST'] = candidates['ST'].fillna('Possible new ST')
    possible_sts = candidates['ST'].to_list()[:9]

    # The original line-by-line writer added an empty trailing field when the
    # row was shorter than 9 columns before appending the possibilities. Keep
    # that here so the final DataFrame slot for "allele_8" stays consistent.
    if len(row_fields) <= 9:
        row_fields.append('')
    row_fields.extend([str(possible_sts), str(possible_alleles)])
    return row_fields


def parse_mlst(mlst_in_path, mlst_out_path, pubmlst_db_root):
    """Parse ``mlst.csv`` row-by-row and write the consolidated TSV.

    The output column count is fixed (see :data:`MLST_HEADER`); short input rows
    are padded and unresolved STs are enriched with candidate alleles from
    PubMLST.
    """
    rows = []
    if os.path.exists(mlst_in_path):
        with open(mlst_in_path, 'r') as fh:
            for line in fh:
                fields = line.strip('\n').strip().split('\t')
                if not fields or not fields[0]:
                    continue
                rows.append(_process_mlst_row(fields, pubmlst_db_root))

    # Pad every row to the header length so the resulting DataFrame is rectangular.
    rows = [r + [''] * (len(MLST_HEADER) - len(r)) for r in rows]
    mlst_df = pd.DataFrame(rows, columns=MLST_HEADER)
    mlst_df.to_csv(mlst_out_path, sep='\t', index=False)


# =============================================================================
# MAIN
# =============================================================================
def main():
    args = get_arguments()
    input_folder = args.input_folder
    out_folder = args.out_dir if args.out_dir else input_folder
    pubmlst_db_root = args.pubMLST_database

    log.info('=== Starting Aluminion Parser ===')
    log.info('Working directory: %s', input_folder)
    preflight_check(input_folder, args)

    # -------------------------------------------------------------------------
    # 1. Independent sub-parsers (phages, integrons, plasmids)
    # -------------------------------------------------------------------------
    if not args.skip_phages:
        log.info('Parsing phages...')
        if hasattr(phage_parser, 'run_parsing'):
            phage_parser.run_parsing(input_folder)

    if not args.skip_integrons:
        log.info('Parsing integrons...')
        if hasattr(integron_parser, 'run_parsing'):
            integron_parser.run_parsing(input_folder)

    if not args.skip_plasmids:
        log.info('Parsing plasmids...')
        if hasattr(copla_parser, 'run_parsing'):
            copla_parser.run_parsing(input_folder, out_folder)

    # -------------------------------------------------------------------------
    # 2. File paths and column schemas
    # -------------------------------------------------------------------------
    mlst_in_path     = os.path.join(input_folder, 'mlst.csv')
    mlst_out_path    = os.path.join(out_folder,   'mlst_modif.csv')
    species_path     = os.path.join(input_folder, '04_taxonomies/kraken2/species.csv')
    genus_path       = os.path.join(input_folder, '04_taxonomies/kraken2/genus.csv')
    gambit_path      = os.path.join(input_folder, '04_taxonomies/gambit.csv')
    kleborate_path   = os.path.join(input_folder, '04_taxonomies/kleborate/enterobacterales__species_output.txt')
    ectyper_path     = os.path.join(input_folder, '04_taxonomies/ectyper/output.tsv')
    abricate_path    = os.path.join(input_folder, 'AbR_report.csv')
    abricate_out     = os.path.join(out_folder,   'AbR_modif.xlsx')
    taxonomy_xlsx    = os.path.join(out_folder,   'taxonomy.xlsx')
    taxonomy_csv     = os.path.join(out_folder,   'taxonomy.csv')

    col_species = ['Sample', 'Perc. reads', 'Species']
    col_genus   = ['Sample', 'Perc. reads', 'Genus']
    kleborate_cols = [
        'strain',
        'enterobacterales__species__species',
        'klebsiella_pneumo_complex__kaptive__K_locus',
        'klebsiella_pneumo_complex__kaptive__K_locus_confidence',
        'klebsiella_pneumo_complex__kaptive__O_locus',
        'klebsiella_pneumo_complex__kaptive__O_locus_confidence',
        'klebsiella_pneumo_complex__amr__Bla_acquired',
        'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired',
        'klebsiella_pneumo_complex__amr__Bla_Carb_acquired',
        'klebsiella_pneumo_complex__virulence_score__virulence_score',
        'klebsiella_pneumo_complex__resistance_score__resistance_score',
        'klebsiella_pneumo_complex__resistance_gene_count__num_resistance_genes',
    ]

    # -------------------------------------------------------------------------
    # 3. Read and reshape per-tool inputs
    # -------------------------------------------------------------------------

    # --- ABRICATE -----------------------------------------------------------
    if not args.skip_abr:
        log.info('Processing resistance genes (Abricate)...')
        abricate_df = safe_read_csv(abricate_path, required_cols=['#FILE'], sep='\t')
        if not abricate_df.empty and '#FILE' in abricate_df.columns:
            # Tag every gene-presence cell with its gene name so the row can be
            # collapsed into a single comma-separated "Resistance_genes" string.
            for col in abricate_df.columns[2:]:
                abricate_df[col] = abricate_df[col].apply(
                    lambda x, c=col: f'{c} ({x})' if x != '.' and not pd.isna(x) else ''
                )

            sample_to_genes = {}
            for i in range(len(abricate_df)):
                row_values = abricate_df.loc[i, :].values.tolist()
                sample_id = row_values[0]
                gene_hits = [v for v in row_values if v != ''][2:]
                sample_to_genes[sample_id] = gene_hits

            abricate_df['Resistance_genes'] = abricate_df['#FILE'].map(sample_to_genes)
            abricate_df['Resistance_genes'] = abricate_df['Resistance_genes'].apply(
                lambda x: ', '.join(x) if isinstance(x, list) else ''
            )
            abricate_df['#FILE'] = (
                abricate_df['#FILE'].astype(str).str.split('/').str[-1].str.rsplit('.', n=1).str[0]
            )
            cols = abricate_df.columns.tolist()
            if len(cols) > 2:
                cols = cols[:2] + cols[-1:] + cols[2:-1]
                abricate_df = abricate_df[cols]
            abricate_df.to_excel(abricate_out, index=False)

    # --- KRAKEN -------------------------------------------------------------
    if not args.skip_kraken:
        log.info('Processing Kraken2 taxonomy...')
        species_df = safe_read_csv(species_path, required_cols=col_species,
                                   sep='\t', names=col_species, header=None, decimal='.')
        genus_df   = safe_read_csv(genus_path,   required_cols=col_genus,
                                   sep='\t', names=col_genus,   header=None, decimal='.')
    else:
        species_df = pd.DataFrame(columns=col_species)
        genus_df   = pd.DataFrame(columns=col_genus)

    # --- TYPING (MLST, Kleborate, GAMBIT, ECTyper) --------------------------
    if not args.skip_typing:
        log.info('Processing typing data (MLST, Kleborate, GAMBIT, ECTyper)...')
        gambit_df    = safe_read_csv(gambit_path,    required_cols=['query', 'closest.description'], sep=',')
        kleborate_df = safe_read_csv(kleborate_path, required_cols=kleborate_cols, sep='\t', usecols=kleborate_cols)
        ectyper_df   = safe_read_csv(ectyper_path,   required_cols=['Name', 'Serotype'], sep='\t', usecols=['Name', 'Serotype'])
        parse_mlst(mlst_in_path, mlst_out_path, pubmlst_db_root)
    else:
        gambit_df    = pd.DataFrame(columns=['query', 'closest.description'])
        kleborate_df = pd.DataFrame(columns=kleborate_cols)
        ectyper_df   = pd.DataFrame(columns=['Name', 'Serotype'])
        # Still emit an empty mlst_modif.csv with the right header so downstream
        # consumers don't crash on a missing file.
        pd.DataFrame(columns=MLST_HEADER).to_csv(mlst_out_path, sep='\t', index=False)

    # -------------------------------------------------------------------------
    # 4. Cross-reference and consolidate
    # -------------------------------------------------------------------------
    log.info('Consolidating data and generating final Excel...')
    mlst_df = safe_read_csv(mlst_out_path, required_cols=['Sample'], sep='\t')
    mlst_df['Sample']      = mlst_df['Sample'].astype(str)
    species_df['Sample']   = species_df['Sample'].astype(str)
    genus_df['Sample']     = genus_df['Sample'].astype(str)
    kleborate_df['strain'] = kleborate_df['strain'].astype(str)
    kleborate_df.fillna('', inplace=True)

    if not gambit_df.empty:
        gambit_df = gambit_df[['query', 'closest.description']]
        gambit_df = gambit_df.rename(columns={'query': 'Sample',
                                              'closest.description': 'Subspecies'})
        # Strip bracketed / parenthesised qualifiers from the GAMBIT subspecies label.
        gambit_df['Subspecies'] = gambit_df['Subspecies'].str.replace(r'\(.*?\)', '', regex=True)
        gambit_df['Subspecies'] = gambit_df['Subspecies'].str.replace(r'\[.*?\]', '', regex=True)
        gambit_df['Sample'] = gambit_df['Sample'].astype(str)
        gambit_df = gambit_df.apply(lambda x: x.str.strip() if x.dtype == 'object' else x)

    if not species_df.empty:
        species_df['Majority_species'] = (
            species_df['Species'] + ' (' +
            species_df['Perc. reads'].fillna(0).round(0).astype(int).astype(str) + '%)'
        )
        # For each sample, keep the first row as the majority species and pack the
        # remaining rows into a comma-separated "Contaminants" string.
        species_df['Contaminants'] = species_df.groupby(['Sample'])['Majority_species'].transform(
            lambda x: '-' if len(x) == 1 else ','.join(x.iloc[1:])
        )
        species_df = species_df.drop_duplicates(subset=['Sample', 'Contaminants'], keep='first')
        species_df.drop(['Perc. reads', 'Species'], axis=1, inplace=True)
    else:
        species_df['Majority_species'] = ''
        species_df['Contaminants']     = ''

    if not genus_df.empty:
        genus_df['Majority_genus'] = (
            genus_df['Genus'] + ' (' +
            genus_df['Perc. reads'].fillna(0).round(0).astype(int).astype(str) + '%)'
        )
        genus_df.drop(['Perc. reads', 'Genus'], axis=1, inplace=True)
    else:
        genus_df['Majority_genus'] = ''

    kraken_df = pd.merge(species_df, genus_df, how='left', on='Sample')
    kraken_cols = ['Sample', 'Majority_genus', 'Majority_species', 'Contaminants']
    for col in kraken_cols:
        if col not in kraken_df.columns:
            kraken_df[col] = ''
    kraken_df = kraken_df[kraken_cols]
    kraken_df['Sample'] = kraken_df['Sample'].astype(str)
    kraken_df = kraken_df.apply(lambda x: x.str.strip() if x.dtype == 'object' else x)

    # Successive left-joins: kraken -> gambit -> mlst -> kleborate -> ectyper
    merged_taxonomy = pd.merge(kraken_df, gambit_df, how='left', on='Sample')
    merged_taxonomy = pd.merge(merged_taxonomy, mlst_df, how='left', on='Sample')

    if not kleborate_df.empty:
        # Strip Kleborate's internal version suffixes (".v1^", "^") and dedupe.
        for amr_col in ('klebsiella_pneumo_complex__amr__Bla_acquired',
                        'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired',
                        'klebsiella_pneumo_complex__amr__Bla_Carb_acquired'):
            kleborate_df[amr_col] = (
                kleborate_df[amr_col]
                .str.replace(r'\.v1\^', '', regex=True)
                .str.replace(r'\^', '', regex=True)
                .str.split(';')
                .apply(lambda x: ', '.join(map(str, set(x))) if isinstance(x, list) else '')
            )

        # Collapse Kaptive "unknown (...)" annotations to a plain "-".
        kleborate_df['klebsiella_pneumo_complex__kaptive__K_locus'] = (
            kleborate_df['klebsiella_pneumo_complex__kaptive__K_locus']
            .str.replace(r'unknown \([A-Z]*[0-9]*\-*[A-Z]*[0-9]*\)', '-', regex=True)
        )
        kleborate_df['klebsiella_pneumo_complex__kaptive__O_locus'] = (
            kleborate_df['klebsiella_pneumo_complex__kaptive__O_locus']
            .str.replace(r'unknown \([A-Z]*[0-9]*\/*[A-Z]*[0-9]*[av]*[0-9]*\)', '-', regex=True)
        )
        kleborate_df['KO_locus'] = (
            kleborate_df['klebsiella_pneumo_complex__kaptive__K_locus'] + '/' +
            kleborate_df['klebsiella_pneumo_complex__kaptive__O_locus']
        )

        kleborate_df = kleborate_df.rename(columns={
            'strain': 'Sample',
            'klebsiella_pneumo_complex__amr__Bla_acquired'                       : 'Other_resistance',
            'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'                  : 'ESBL',
            'klebsiella_pneumo_complex__amr__Bla_Carb_acquired'                  : 'Carbapenemase',
            'klebsiella_pneumo_complex__virulence_score__virulence_score'        : 'VIRscore',
            'klebsiella_pneumo_complex__resistance_score__resistance_score'      : 'AMRscore',
            'klebsiella_pneumo_complex__resistance_gene_count__num_resistance_genes': 'N_AMR_genes',
        })
        kleborate_df = kleborate_df[['Sample', 'KO_locus', 'Carbapenemase', 'ESBL',
                                     'Other_resistance', 'N_AMR_genes', 'AMRscore', 'VIRscore']]
    else:
        kleborate_df = pd.DataFrame(columns=['Sample', 'KO_locus', 'Carbapenemase', 'ESBL',
                                             'Other_resistance', 'N_AMR_genes', 'AMRscore',
                                             'VIRscore'])

    merged_taxonomy = pd.merge(merged_taxonomy, kleborate_df, how='left', on='Sample')

    ectyper_df = ectyper_df.rename(columns={'Name': 'Sample'})
    ectyper_df['Sample'] = ectyper_df['Sample'].astype(str)
    final_taxonomy = pd.merge(merged_taxonomy, ectyper_df, how='left', on='Sample')

    final_columns = [
        'Sample', 'Majority_genus', 'Majority_species', 'Subspecies',
        'MLST', 'Serotype', 'KO_locus', 'Contaminants',
        'Carbapenemase', 'ESBL', 'Other_resistance',
        'N_AMR_genes', 'AMRscore', 'VIRscore',
        'MLST_scheme', 'allele_1', 'allele_2', 'allele_3', 'allele_4',
        'allele_5', 'allele_6', 'allele_7',
        'Possible_MLSTs', 'Possible_alleles',
    ]
    for col in final_columns:
        if col not in final_taxonomy.columns:
            final_taxonomy[col] = ''

    final_taxonomy = final_taxonomy[final_columns]
    final_taxonomy.to_excel(taxonomy_xlsx, index=False)
    final_taxonomy.to_csv(taxonomy_csv, index=False)

    # -------------------------------------------------------------------------
    # 5. Kraken + MLST quick-reference (kraken_mlst.xlsx)
    # -------------------------------------------------------------------------
    kraken_mlst_df = final_taxonomy[[
        'Sample', 'Majority_genus', 'Majority_species', 'Subspecies',
        'MLST_scheme', 'MLST', 'allele_1', 'allele_2', 'allele_3', 'allele_4',
        'allele_5', 'allele_6', 'allele_7', 'Possible_MLSTs', 'Possible_alleles',
    ]].copy()

    if not kraken_mlst_df.empty:
        # Split "Genus (NN%)" into its two visible parts for the quick-ref view.
        genus_split = kraken_mlst_df['Majority_genus'].astype(str).str.split(' ', n=1, expand=True)
        kraken_mlst_df['Genus']     = genus_split[0] if 0 in genus_split.columns else ''
        kraken_mlst_df['Genus_pct'] = (
            genus_split[1].str.strip().str.rstrip(')').str.lstrip('(')
            if 1 in genus_split.columns else ''
        )
        species_split = kraken_mlst_df['Majority_species'].astype(str).str.rsplit(' ', n=1, expand=True)
        kraken_mlst_df['Species']     = species_split[0] if 0 in species_split.columns else ''
        kraken_mlst_df['Species_pct'] = (
            species_split[1].str.strip().str.rstrip(')').str.lstrip('(')
            if 1 in species_split.columns else ''
        )
    else:
        kraken_mlst_df['Genus']       = ''
        kraken_mlst_df['Genus_pct']   = ''
        kraken_mlst_df['Species']     = ''
        kraken_mlst_df['Species_pct'] = ''

    kraken_mlst_df.drop(columns=['Majority_genus', 'Majority_species'], inplace=True, errors='ignore')
    kraken_mlst_df['Contaminant'] = ''
    kraken_mlst_df = kraken_mlst_df[[
        'Sample', 'Genus_pct', 'Genus', 'Species_pct', 'Species', 'Contaminant',
        'Subspecies', 'MLST_scheme', 'MLST',
        'allele_1', 'allele_2', 'allele_3', 'allele_4',
        'allele_5', 'allele_6', 'allele_7',
        'Possible_MLSTs', 'Possible_alleles',
    ]]

    kraken_mlst_df.to_excel(os.path.join(out_folder, 'kraken_mlst.xlsx'), index=False)

    log.info('Process completed successfully.')


if __name__ == '__main__':
    main()
