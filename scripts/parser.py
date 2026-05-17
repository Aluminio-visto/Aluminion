#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#############################################################
# Jorge R Grande - HUMV - IDIVAL - Santander
# parser.py takes in a folder containing the output of aluminion.sh
# and outputs a short report summarizing stats from the MinION run
#############################################################

import os
import pandas as pd
import numpy as np
import argparse
import warnings

# ==========================================
# IMPORT EXTERNAL MODULES
# ==========================================
try:
    import phage_parser
    import integron_parser
    import copla_parser
    import Datos_seq_unified2 as run_info_parser
except ImportError as e:
    warnings.warn(f"Missing script in scripts/ folder: {e}")

# ==========================================
# HELPER FUNCTIONS
# ==========================================
def safe_read_csv(filepath, required_cols, sep=',', **kwargs):
    """Safely reads a CSV, returning empty DataFrame with required columns on failure."""
    if not os.path.exists(filepath):
        print(f"\033[93m[WARNING]\033[0m File not found, skipping: {filepath}")
        return pd.DataFrame(columns=required_cols)
    try:
        if os.stat(filepath).st_size == 0:
            print(f"\033[93m[WARNING]\033[0m Empty file: {filepath}")
            return pd.DataFrame(columns=required_cols)

        df = pd.read_csv(filepath, sep=sep, **kwargs)
        for col in required_cols:
            if col not in df.columns:
                df[col] = pd.Series(dtype='object')
        return df
    except Exception as e:
        print(f"\033[93m[WARNING]\033[0m Error reading {filepath}: {e}")
        return pd.DataFrame(columns=required_cols)

def get_arguments():
    parser = argparse.ArgumentParser(prog='parser.py', description='Aluminion Parser: Generates aggregated WGS reports.')
    
    input_group = parser.add_argument_group('Input', 'Input parameters')
    input_group.add_argument('-i', '--input_folder', dest="input_folder", required=True, help="Required. Input folder to be summarized", type=os.path.abspath)
    input_group.add_argument('-db', '--pubMLST_database', dest="pubMLST_database", required=False, default="/home/usuario/miniconda3/envs/mlst/db/pubmlst/", help="Folder containing MLST schemas", type=os.path.abspath)
    
    output_group = parser.add_argument_group('Output', 'Output parameters')
    output_group.add_argument('-o', '--out_dir', dest='out_dir', required=False, help='Final output folder', type=os.path.abspath)
    
    # SKIP FLAGS
    skip_group = parser.add_argument_group('Skip Modules', 'Flags to skip specific analysis parsing')
    skip_group.add_argument('--skip-phages', action='store_true', help="Skip phage parsing (Phastest)")
    skip_group.add_argument('--skip-integrons', action='store_true', help="Skip integron parsing (Integron Finder)")
    skip_group.add_argument('--skip-plasmids', action='store_true', help="Skip plasmid parsing (Copla)")
    skip_group.add_argument('--skip-typing', action='store_true', help="Skip typing data (MLST, Kleborate, Ectyper, Gambit)")
    skip_group.add_argument('--skip-kraken', action='store_true', help="Skip Kraken2 contamination parsing")
    skip_group.add_argument('--skip-abr', action='store_true', help="Skip resistance gene extraction (Abricate)")

    # EXTRA METADATA
    parser.add_argument('--include-run-info', type=str, default=None, help="Path to MinKNOW final_summary.txt file")
    
    return parser.parse_args()


def preflight_check(input_folder, args):
    """Warn about missing input files before processing begins."""
    checks = []

    if not args.skip_abr:
        checks.append(("AbR_report.csv",        os.path.join(input_folder, "AbR_report.csv"),        "--skip-abr"))
    if not args.skip_kraken:
        checks.append(("kraken2/species.csv",    os.path.join(input_folder, "04_taxonomies/kraken2/species.csv"), "--skip-kraken"))
        checks.append(("kraken2/genus.csv",      os.path.join(input_folder, "04_taxonomies/kraken2/genus.csv"),   "--skip-kraken"))
    if not args.skip_typing:
        checks.append(("gambit.csv",             os.path.join(input_folder, "04_taxonomies/gambit.csv"),          "--skip-typing"))
        checks.append(("kleborate output",       os.path.join(input_folder, "04_taxonomies/kleborate/enterobacterales__species_output.txt"), "--skip-typing"))
        checks.append(("ectyper output.tsv",     os.path.join(input_folder, "04_taxonomies/ectyper/output.tsv"),  "--skip-typing"))
        checks.append(("mlst.csv",               os.path.join(input_folder, "mlst.csv"),                          "--skip-typing"))

    missing = [(name, flag) for name, path, flag in checks if not os.path.exists(path)]
    if missing:
        print("\033[93m[PREFLIGHT]\033[0m The following input files were not found:")
        for name, flag in missing:
            print(f"  \033[93m✗\033[0m  {name}  (skip with {flag})")
        print("\033[93m[PREFLIGHT]\033[0m Processing will continue — missing modules will produce empty tables.\033[0m\n")
    else:
        print("\033[92m[PREFLIGHT]\033[0m All expected input files found.\033[0m\n")


def main():
    args = get_arguments()
    input_folder = args.input_folder
    out_folder = args.out_dir if args.out_dir else input_folder
    DB = args.pubMLST_database

    print(f"\n\033[94m=== Starting Aluminion Parser ===\033[0m")
    print(f"Working directory: {input_folder}\n")
    preflight_check(input_folder, args)

    # ==========================================
    # 1. RUN INDEPENDENT MODULES
    # ==========================================
    if not args.skip_phages:
        print("-> Parsing phages...")
        if hasattr(phage_parser, 'run_parsing'): phage_parser.run_parsing(input_folder)

    if not args.skip_integrons:
        print("-> Parsing integrons...")
        if hasattr(integron_parser, 'run_parsing'): integron_parser.run_parsing(input_folder)

    if not args.skip_plasmids:
        print("-> Parsing plasmids...")
        if hasattr(copla_parser, 'run_parsing'): copla_parser.run_parsing(input_folder, out_folder)

    # ==========================================
    # 2. PATH AND COLUMN DEFINITIONS
    # ==========================================
    report              = input_folder + "/mlst.csv"
    report_2            = out_folder   + "/mlst_modif.csv"
    species_report      = input_folder + "/04_taxonomies/kraken2/species.csv"
    genus_report        = input_folder + "/04_taxonomies/kraken2/genus.csv"
    gambit_report       = input_folder + "/04_taxonomies/gambit.csv"
    kleb_report         = input_folder + "/04_taxonomies/kleborate/enterobacterales__species_output.txt"
    ec                  = input_folder + "/04_taxonomies/ectyper/output.tsv"
    abricate            = input_folder + "/AbR_report.csv"
    abricate_out        = out_folder   + "/AbR_modif.xlsx"
    taxonomy_out        = out_folder   + "/taxonomy.xlsx"
    taxonomy_csv        = out_folder   + "/taxonomy.csv"

    col_species = ['Sample','Perc. reads', 'Species']
    col_generos = ['Sample','Perc. reads', 'Genus']
    kleb_cols = ['strain','enterobacterales__species__species','klebsiella_pneumo_complex__kaptive__K_locus',
                 'klebsiella_pneumo_complex__kaptive__K_locus_confidence','klebsiella_pneumo_complex__kaptive__O_locus',
                 'klebsiella_pneumo_complex__kaptive__O_locus_confidence','klebsiella_pneumo_complex__amr__Bla_acquired',
                 'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired','klebsiella_pneumo_complex__amr__Bla_Carb_acquired',
                 'klebsiella_pneumo_complex__virulence_score__virulence_score',
                 'klebsiella_pneumo_complex__resistance_score__resistance_score',
                 'klebsiella_pneumo_complex__resistance_gene_count__num_resistance_genes']

    # ==========================================
    # 3. READ INPUT FILES
    # ==========================================

    # ABRICATE
    if not args.skip_abr:
        print("-> Processing resistance genes (Abricate)...")
        df_abricate = safe_read_csv(abricate, required_cols=['#FILE'], sep='\t')
        if not df_abricate.empty and '#FILE' in df_abricate.columns:
            for col in df_abricate.columns[2:]:
                df_abricate[col] = df_abricate[col].apply(lambda x: col + ' (' + str(x) +')' if x != '.' and not pd.isna(x) else '')
            diccionario={}
            for i in range(len(df_abricate)):
                value = (df_abricate.loc[i, :].values.tolist())
                muestra = value[0]
                valores = [valor for valor in value if valor != ''][2:]
                diccionario[muestra] = valores

            df_abricate["Resistance_genes"] = df_abricate["#FILE"].map(diccionario)
            df_abricate["Resistance_genes"] = df_abricate["Resistance_genes"].apply(lambda x: ', '.join([i for i in x]) if isinstance(x, list) else "")
            df_abricate["#FILE"] = df_abricate["#FILE"].astype(str).str.split('/').str[-1].str.rsplit('.', n=1).str[0]
            cols = df_abricate.columns.tolist()
            if len(cols) > 2:
                cols = cols[:2] + cols[-1:] + cols[2:-1]
                df_abricate  = df_abricate[cols]
            df_abricate.to_excel(abricate_out, index=False)

    # KRAKEN
    if not args.skip_kraken:
        print("-> Processing Kraken2 taxonomy...")
        species = safe_read_csv(species_report, required_cols=col_species, sep='\t', names=col_species, header=None, decimal='.')
        genus   = safe_read_csv(genus_report, required_cols=col_generos, sep='\t', names=col_generos, header=None, decimal='.')
    else:
        species = pd.DataFrame(columns=col_species)
        genus   = pd.DataFrame(columns=col_generos)

    # TYPING (MLST, Kleborate, Gambit, Ectyper)
    if not args.skip_typing:
        print("-> Processing typing data (MLST, Kleborate, Gambit, Ectyper)...")
        gambit   = safe_read_csv(gambit_report, required_cols=['query', 'closest.description'], sep=',')
        kleborate = safe_read_csv(kleb_report, required_cols=kleb_cols, sep='\t', usecols=kleb_cols)
        ectyper  = safe_read_csv(ec, required_cols=['Name','Serotype'], sep='\t', usecols=['Name','Serotype'])
        
        # MLST processing
        cabecera = ['Sample','MLST_scheme','MLST','allele_1','allele_2','allele_3','allele_4','allele_5','allele_6','allele_7','allele_8','Possible_MLSTs','Possible_alleles']
        with open(report_2, 'w') as f:
            f.write('\t'.join(map(str, cabecera)) +'\n')

        if os.path.exists(report):
            with open(report, 'r') as file:
                for line in file:
                    column = line.strip('\n').strip().split('\t')
                    organismo = column[1]
                    muestra   = column[0].split('/')[1].split('.')[0] if '/' in column[0] else column[0]
                    column[0] = muestra
                    db_org    = DB + '/' + organismo + '/' + organismo + '.txt'

                    if column[1] == '-':
                        column[1] = 'No associated MLST scheme'
                        with open(report_2, 'a') as f:
                            f.write('\t'.join(map(str, column)) +'\n')            
                    elif column[1] != '-' and column[2] != '-':
                        with open(report_2, 'a') as f:
                            f.write('\t'.join(map(str, column)) +'\n')          
                    elif column[1] != '-' and column[2] == '-':
                        genes   = [col.split('(')[0] for col in column[3:]]
                        numeros = [col.split('(')[1].strip(')') for col in column[3:]]
                        diccion = dict(zip(genes, numeros))
                        diccion = {k: v for k, v in diccion.items() if v.isdigit()}
                       
                        mlst_df = pd.DataFrame([diccion])
                        mlst_df = mlst_df.astype('int64')

                        if os.path.exists(db_org):
                            df = pd.read_csv(db_org, sep = '\t')
                        else:
                            df = pd.DataFrame(columns=list(diccion.keys()) + ['ST'])

                        genes_seguros = list(diccion.keys())
                        if len(genes_seguros) == 0:
                            posibles = pd.DataFrame([{'ST': []}])
                        else:
                            posibles = pd.merge(mlst_df, df, how='left', on=genes_seguros)

                        if posibles.shape[0] == 1 and posibles['ST'].isna().any():
                            posibles_alelos = posibles.iloc[:,posibles.columns.get_indexer(['ST'])]
                        else:
                            posibles_alelos = posibles.iloc[:,posibles.columns.get_indexer(['ST'])+1]

                        posibles_alelos = posibles_alelos.to_dict('list')
                        for k, v in posibles_alelos.items():
                            v = ['No allele for this gene matches this combination' if pd.isna(x) else x for x in v]
                            posibles_alelos[k] = v[:9]

                        posibles['ST'].fillna('Possible new ST', inplace=True)
                        posibles_ST = posibles.ST.to_list()[:9]
                        if len(column) <= 9:
                            column.append(str(''))     
                        column.append(posibles_ST)
                        column.append(posibles_alelos)

                        with open(report_2, 'a') as f:
                            f.write('\t'.join(map(str, column)) +'\n')
    else:
        gambit = pd.DataFrame(columns=['query', 'closest.description'])
        kleborate = pd.DataFrame(columns=kleb_cols)
        ectyper = pd.DataFrame(columns=['Name','Serotype'])
        # Create empty mlst_modif
        cabecera = ['Sample','MLST_scheme','MLST','allele_1','allele_2','allele_3','allele_4','allele_5','allele_6','allele_7','allele_8','Possible_MLSTs','Possible_alleles']
        with open(report_2, 'w') as f:
            f.write('\t'.join(map(str, cabecera)) +'\n')


    # ==========================================
    # 4. CROSS-REFERENCE AND CONSOLIDATION (PANDAS MERGE)
    # ==========================================
    print("-> Consolidating data and generating final Excel...")
    mlst = safe_read_csv(report_2, required_cols=['Sample'], sep='\t')
    mlst['Sample'] = mlst['Sample'].astype(str)
    species['Sample'] = species['Sample'].astype(str)
    genus['Sample']   = genus['Sample'].astype(str)
    kleborate['strain'] = kleborate['strain'].astype(str)
    kleborate.fillna('', inplace=True)

    if not gambit.empty:
        gambit = gambit[['query','closest.description']]
        gambit = gambit.rename(columns={'query':'Sample','closest.description':'Subspecies'})
        gambit['Subspecies'] = gambit['Subspecies'].str.replace(r"\(.*?\)", "", regex=True)
        gambit['Subspecies'] = gambit['Subspecies'].str.replace(r"\[.*?\]", "", regex=True)
        gambit['Sample'] = gambit['Sample'].astype(str)
        gambit = gambit.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    if not species.empty:
        species['Majority_species'] = species['Species'] + ' (' + species["Perc. reads"].fillna(0).round(0).astype(int).astype(str) +  '%)'
        species['Contaminants'] = species.groupby(['Sample'])['Majority_species'].transform(lambda x: '-' if len(x) == 1 else ','.join(x.iloc[1:]))
        species = species.drop_duplicates(subset=['Sample','Contaminants'], keep = 'first')
        species.drop(['Perc. reads','Species'], axis = 1, inplace=True)
    else:
        species['Majority_species'] = ""
        species['Contaminants'] = ""

    if not genus.empty:
        genus['Majority_genus']  = genus['Genus'] + ' (' + genus['Perc. reads'].fillna(0).round(0).astype(int).astype(str) +  '%)'
        genus.drop(['Perc. reads','Genus'], axis = 1, inplace=True)
    else:
        genus['Majority_genus'] = ""

    kraken = pd.merge(species, genus, how = "left", on='Sample')
    kraken_cols = ['Sample','Majority_genus','Majority_species','Contaminants']
    for col in kraken_cols:
        if col not in kraken.columns: kraken[col] = ""
    kraken = kraken[kraken_cols]
    kraken['Sample'] = kraken['Sample'].astype(str)
    kraken = kraken.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    subesp = pd.merge(kraken, gambit, how="left", on='Sample')
    intermedio = pd.merge(subesp, mlst, how='left', on='Sample')

    if not kleborate.empty:
        kleborate['klebsiella_pneumo_complex__amr__Bla_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_acquired'].str.replace(r'\.v1\^', '', regex=True).str.replace(r'\^', '', regex=True).str.split(';').apply(lambda x: ', '.join(map(str, set(x))) if isinstance(x, list) else "")
        kleborate['klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'].str.replace(r'\.v1\^', '', regex=True).str.replace(r'\^', '', regex=True).str.split(';').apply(lambda x: ', '.join(map(str, set(x))) if isinstance(x, list) else "")
        kleborate['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'].str.replace(r'\.v1\^', '', regex=True).str.replace(r'\^', '', regex=True).str.split(';').apply(lambda x: ', '.join(map(str, set(x))) if isinstance(x, list) else "")
        
        kleborate['klebsiella_pneumo_complex__kaptive__K_locus'] = kleborate['klebsiella_pneumo_complex__kaptive__K_locus'].str.replace(r'unknown \([A-Z]*[0-9]*\-*[A-Z]*[0-9]*\)', '-', regex=True)  
        kleborate['klebsiella_pneumo_complex__kaptive__O_locus'] = kleborate['klebsiella_pneumo_complex__kaptive__O_locus'].str.replace(r'unknown \([A-Z]*[0-9]*\/*[A-Z]*[0-9]*[av]*[0-9]*\)', '-', regex=True)
        kleborate['KO_locus'] = kleborate['klebsiella_pneumo_complex__kaptive__K_locus'] + '/' + kleborate['klebsiella_pneumo_complex__kaptive__O_locus']

        kleborate = kleborate.rename(columns={'strain':'Sample',
                            'klebsiella_pneumo_complex__amr__Bla_acquired':'Other_resistance',
                            'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired':'ESBL',
                            'klebsiella_pneumo_complex__amr__Bla_Carb_acquired':'Carbapenemase',
                            'klebsiella_pneumo_complex__virulence_score__virulence_score':'VIRscore',
                            'klebsiella_pneumo_complex__resistance_score__resistance_score':'AMRscore',
                            'klebsiella_pneumo_complex__resistance_gene_count__num_resistance_genes':'N_AMR_genes'})
        kleborate = kleborate[['Sample', 'KO_locus','Carbapenemase','ESBL','Other_resistance', 'N_AMR_genes','AMRscore','VIRscore']]
    else:
        kleborate = pd.DataFrame(columns=['Sample', 'KO_locus','Carbapenemase','ESBL','Other_resistance', 'N_AMR_genes','AMRscore','VIRscore'])

    resultado = pd.merge(intermedio, kleborate, how="left", on='Sample')

    ectyper = ectyper.rename(columns={'Name':'Sample'})
    ectyper['Sample'] = ectyper['Sample'].astype(str)
    resultado_final = pd.merge(resultado, ectyper, how="left", on='Sample')

    # Optional MinKNOW metadata integration
    if args.include_run_info and os.path.exists(args.include_run_info):
        print(f"-> Integrating MinKNOW metadata from: {args.include_run_info}")
        if hasattr(run_info_parser, 'parse_minion_sum'):
            
            pass

    columnas_finales = ['Sample','Majority_genus','Majority_species','Subspecies','MLST','Serotype','KO_locus','Contaminants','Carbapenemase','ESBL','Other_resistance','N_AMR_genes','AMRscore','VIRscore','MLST_scheme','allele_1','allele_2','allele_3','allele_4','allele_5','allele_6','allele_7','Possible_MLSTs','Possible_alleles']
    for col in columnas_finales:
        if col not in resultado_final.columns: resultado_final[col] = ""

    resultado_final = resultado_final[columnas_finales]
    resultado_final.to_excel(taxonomy_out, index=False)
    resultado_final.to_csv(taxonomy_csv, index=False)

    ultimo_df = resultado_final[['Sample','Majority_genus','Majority_species','Subspecies','MLST_scheme','MLST','allele_1','allele_2','allele_3','allele_4','allele_5','allele_6','allele_7','Possible_MLSTs','Possible_alleles']].copy()

    if not ultimo_df.empty:
        provis1 = ultimo_df['Majority_genus'].astype(str).str.split(" ", n=1, expand=True)
        ultimo_df["Genus"]     = provis1[0] if 0 in provis1.columns else ""
        ultimo_df["Genus_pct"] = provis1[1].str.strip().str.rstrip(')').str.lstrip('(') if 1 in provis1.columns else ""
        provis2 = ultimo_df['Majority_species'].astype(str).str.rsplit(" ", n=1, expand=True)
        ultimo_df["Species"]    = provis2[0] if 0 in provis2.columns else ""
        ultimo_df["Species_pct"] = provis2[1].str.strip().str.rstrip(')').str.lstrip('(') if 1 in provis2.columns else ""
    else:
        ultimo_df["Genus"] = ""; ultimo_df["Genus_pct"] = ""; ultimo_df["Species"] = ""; ultimo_df["Species_pct"] = ""

    ultimo_df.drop(columns=['Majority_genus', 'Majority_species'], inplace=True, errors='ignore')
    ultimo_df["Contaminant"] = ""
    ultimo_df = ultimo_df[['Sample',"Genus_pct", "Genus","Species_pct", "Species", "Contaminant",'Subspecies','MLST_scheme','MLST','allele_1','allele_2','allele_3','allele_4','allele_5','allele_6','allele_7','Possible_MLSTs','Possible_alleles']]

    ultimo_df.to_excel(os.path.join(out_folder, "kraken_mlst.xlsx"), index=False)

    print("\n\033[92m[INFO]\033[0m Process completed successfully!")

if __name__ == '__main__':
    main()