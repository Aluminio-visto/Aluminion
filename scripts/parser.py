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
# IMPORTACIÓN DE MÓDULOS EXTERNOS
# ==========================================
try:
    import phage_parser
    import integron_parser
    import copla_parser
    import Datos_seq_unified2 as run_info_parser
except ImportError as e:
    warnings.warn(f"Falta algún script en la carpeta scripts/: {e}")

# ==========================================
# FUNCIONES AUXILIARES
# ==========================================
def safe_read_csv(filepath, required_cols, sep=',', **kwargs):
    """Lee un CSV de forma segura devolviendo columnas vacías si falla."""
    if not os.path.exists(filepath):
        print(f"\033[93m[WARNING]\033[0m Archivo no encontrado. Saltando: {filepath}")
        return pd.DataFrame(columns=required_cols)
    try:
        if os.stat(filepath).st_size == 0:
            print(f"\033[93m[WARNING]\033[0m Archivo vacío: {filepath}")
            return pd.DataFrame(columns=required_cols)
            
        df = pd.read_csv(filepath, sep=sep, **kwargs)
        for col in required_cols:
            if col not in df.columns:
                df[col] = pd.Series(dtype='object')
        return df
    except Exception as e:
        print(f"\033[93m[WARNING]\033[0m Error al leer {filepath}: {e}")
        return pd.DataFrame(columns=required_cols)

def get_arguments():
    parser = argparse.ArgumentParser(prog='parser.py', description='Aluminion Parser: Generates aggregated WGS reports.')
    
    input_group = parser.add_argument_group('Input', 'Input parameters')
    input_group.add_argument('-i', '--input_folder', dest="input_folder", required=True, help="Required. Input folder to be summarized", type=os.path.abspath)
    input_group.add_argument('-db', '--pubMLST_database', dest="pubMLST_database", required=False, default="/home/usuario/miniconda3/envs/mlst/db/pubmlst/", help="Folder containing MLST schemas", type=os.path.abspath)
    
    output_group = parser.add_argument_group('Output', 'Output parameters')
    output_group.add_argument('-o', '--out_dir', dest='out_dir', required=False, help='Final output folder', type=os.path.abspath)
    
    # NUEVOS PARÁMETROS DE OMISIÓN
    skip_group = parser.add_argument_group('Skip Modules', 'Flags to skip specific analysis parsing')
    skip_group.add_argument('--skip-phages', action='store_true', help="Omite el parseo de fagos (Phastest)")
    skip_group.add_argument('--skip-integrons', action='store_true', help="Omite el parseo de integrones (Integron Finder)")
    skip_group.add_argument('--skip-plasmids', action='store_true', help="Omite la lectura de plásmidos (Copla)")
    skip_group.add_argument('--skip-typing', action='store_true', help="Omite la lectura de tipado (MLST, Kleborate, Ectyper, Gambit)")
    skip_group.add_argument('--skip-kraken', action='store_true', help="Omite la lectura de contaminaciones de Kraken2")
    skip_group.add_argument('--skip-abr', action='store_true', help="Omite la extracción de genes de resistencia (Abricate)")
    
    # METADATOS EXTRA
    parser.add_argument('--include-run-info', type=str, default=None, help="Ruta al archivo final_summary.txt de MinKNOW")
    
    return parser.parse_args()


def main():
    args = get_arguments()
    input_folder = args.input_folder
    out_folder = args.out_dir if args.out_dir else input_folder
    DB = args.pubMLST_database

    print(f"\n\033[94m=== Iniciando Aluminion Parser ===\033[0m")
    print(f"Directorio de trabajo: {input_folder}\n")

    # ==========================================
    # 1. EJECUCIÓN DE MÓDULOS INDEPENDIENTES
    # ==========================================
    if not args.skip_phages:
        print("-> Parseando Fagos...")
        if hasattr(phage_parser, 'run_parsing'): phage_parser.run_parsing(input_folder)
            
    if not args.skip_integrons:
        print("-> Parseando Integrones...")
        if hasattr(integron_parser, 'run_parsing'): integron_parser.run_parsing(input_folder)

    if not args.skip_plasmids:
        print("-> Parseando Plásmidos...")
        if hasattr(copla_parser, 'run_parsing'): copla_parser.run_parsing(input_folder, out_folder)

    # ==========================================
    # 2. DEFINICIÓN DE RUTAS Y COLUMNAS
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

    col_species = ['Sample','Perc. reads', 'Especie']
    col_generos = ['Sample','Perc. reads', 'Género']
    kleb_cols = ['strain','enterobacterales__species__species','klebsiella_pneumo_complex__kaptive__K_locus',
                 'klebsiella_pneumo_complex__kaptive__K_locus_confidence','klebsiella_pneumo_complex__kaptive__O_locus',
                 'klebsiella_pneumo_complex__kaptive__O_locus_confidence','klebsiella_pneumo_complex__amr__Bla_acquired',
                 'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired','klebsiella_pneumo_complex__amr__Bla_Carb_acquired',
                 'klebsiella_pneumo_complex__virulence_score__virulence_score',
                 'klebsiella_pneumo_complex__resistance_score__resistance_score',
                 'klebsiella_pneumo_complex__resistance_gene_count__num_resistance_genes']

    # ==========================================
    # 3. LECTURA DE ARCHIVOS
    # ==========================================
    
    # ABRICATE
    if not args.skip_abr:
        print("-> Procesando genes de resistencia (Abricate)...")
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

            df_abricate["Genes resistencia"] = df_abricate["#FILE"].map(diccionario)
            df_abricate["Genes resistencia"] = df_abricate["Genes resistencia"].apply(lambda x: ', '.join([i for i in x]) if isinstance(x, list) else "")
            df_abricate["#FILE"] = df_abricate["#FILE"].astype(str).str.split('/').str[-1]
            cols = df_abricate.columns.tolist()
            if len(cols) > 2:
                cols = cols[:2] + cols[-1:] + cols[2:-1]
                df_abricate  = df_abricate[cols]
            df_abricate.to_excel(abricate_out, index=False)

    # KRAKEN
    if not args.skip_kraken:
        print("-> Procesando taxonomía de Kraken2...")
        species = safe_read_csv(species_report, required_cols=col_species, sep='\t', names=col_species, header=None, decimal='.')
        genus   = safe_read_csv(genus_report, required_cols=col_generos, sep='\t', names=col_generos, header=None, decimal='.')
    else:
        species = pd.DataFrame(columns=col_species)
        genus   = pd.DataFrame(columns=col_generos)

    # TIPADO (MLST, Kleborate, Gambit, Ectyper)
    if not args.skip_typing:
        print("-> Procesando datos de tipado (MLST, Kleborate, Gambit, Ectyper)...")
        gambit   = safe_read_csv(gambit_report, required_cols=['query', 'closest.description'], sep=',')
        kleborate = safe_read_csv(kleb_report, required_cols=kleb_cols, sep='\t', usecols=kleb_cols)
        ectyper  = safe_read_csv(ec, required_cols=['Name','Serotype'], sep='\t', usecols=['Name','Serotype'])
        
        # Procesamiento MLST
        cabecera = ['Sample','Esquema MLST','MLST','alelo #1','alelo #2','alelo #3','alelo #4','alelo #5','alelo #6','alelo #7','alelo #8','MLSTs posibles','Alelos posibles']
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
                        column[1] = 'No tiene esquema asociado'
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
                            v = ['Ningún alelo de este gen coincide con esta combinación' if pd.isna(x) else x for x in v]
                            posibles_alelos[k] = v[:9]                    

                        posibles['ST'].fillna('Posible nuevo ST', inplace=True)
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
        # Creamos el mlst_modif vacío
        cabecera = ['Sample','Esquema MLST','MLST','alelo #1','alelo #2','alelo #3','alelo #4','alelo #5','alelo #6','alelo #7','alelo #8','MLSTs posibles','Alelos posibles']
        with open(report_2, 'w') as f:
            f.write('\t'.join(map(str, cabecera)) +'\n')


    # ==========================================
    # 4. CRUCE Y CONSOLIDACIÓN (MERGE DE PANDAS)
    # ==========================================
    print("-> Consolidando datos y generando excel final...")
    mlst = safe_read_csv(report_2, required_cols=['Sample'], sep='\t')
    mlst['Sample'] = mlst['Sample'].astype(str)
    species['Sample'] = species['Sample'].astype(str)
    genus['Sample']   = genus['Sample'].astype(str)
    kleborate['strain'] = kleborate['strain'].astype(str)
    kleborate.fillna('', inplace=True)

    if not gambit.empty:
        gambit = gambit[['query','closest.description']]
        gambit = gambit.rename(columns={'query':'Sample','closest.description':'Subespecie'})
        gambit['Subespecie'] = gambit['Subespecie'].str.replace(r"\(.*?\)", "", regex=True)
        gambit['Subespecie'] = gambit['Subespecie'].str.replace(r"\[.*?\]", "", regex=True)
        gambit['Sample'] = gambit['Sample'].astype(str)
        gambit = gambit.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    if not species.empty:
        species['Especie mayoritaria'] = species['Especie'] + ' (' + species["Perc. reads"].fillna(0).round(0).astype(int).astype(str) +  '%)'
        species['Posibles contaminantes'] = species.groupby(['Sample'])['Especie mayoritaria'].transform(lambda x: '-' if len(x) == 1 else ','.join(x.iloc[1:]))
        species = species.drop_duplicates(subset=['Sample','Posibles contaminantes'], keep = 'first')
        species.drop(['Perc. reads','Especie'], axis = 1, inplace=True)
    else:
        species['Especie mayoritaria'] = ""
        species['Posibles contaminantes'] = ""

    if not genus.empty:
        genus['Género mayoritario']  = genus['Género'] + ' (' + genus['Perc. reads'].fillna(0).round(0).astype(int).astype(str) +  '%)'
        genus.drop(['Perc. reads','Género'], axis = 1, inplace=True)
    else:
        genus['Género mayoritario'] = ""

    kraken = pd.merge(species, genus, how = "left", on='Sample')
    kraken_cols = ['Sample','Género mayoritario','Especie mayoritaria','Posibles contaminantes']
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
        kleborate['K/O locus'] = kleborate['klebsiella_pneumo_complex__kaptive__K_locus'] + '/' + kleborate['klebsiella_pneumo_complex__kaptive__O_locus']
        
        kleborate = kleborate.rename(columns={'strain':'Sample',
                            'klebsiella_pneumo_complex__amr__Bla_acquired':'Otras',
                            'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired':'BLEE adquirida',
                            'klebsiella_pneumo_complex__amr__Bla_Carb_acquired':'Carba adquirida',
                            'klebsiella_pneumo_complex__virulence_score__virulence_score':'VIRscore',
                            'klebsiella_pneumo_complex__resistance_score__resistance_score':'AMRscore',
                            'klebsiella_pneumo_complex__resistance_gene_count__num_resistance_genes':'Nº genes AMR'})
        kleborate = kleborate[['Sample', 'K/O locus','Carba adquirida','BLEE adquirida','Otras', 'Nº genes AMR','AMRscore','VIRscore']]
    else:
        kleborate = pd.DataFrame(columns=['Sample', 'K/O locus','Carba adquirida','BLEE adquirida','Otras', 'Nº genes AMR','AMRscore','VIRscore'])

    resultado = pd.merge(intermedio, kleborate, how="left", on='Sample')

    ectyper = ectyper.rename(columns={'Name':'Sample'})
    ectyper['Sample'] = ectyper['Sample'].astype(str)
    resultado_final = pd.merge(resultado, ectyper, how="left", on='Sample')

    # Integración info MinKNOW opcional
    if args.include_run_info and os.path.exists(args.include_run_info):
        print(f"-> Integrando metadatos de MinKNOW desde: {args.include_run_info}")
        if hasattr(run_info_parser, 'parse_minion_sum'):
            
            pass

    columnas_finales = ['Sample','Género mayoritario','Especie mayoritaria','Subespecie','MLST','Serotype','K/O locus','Posibles contaminantes','Carba adquirida','BLEE adquirida','Otras','Nº genes AMR','AMRscore','VIRscore','Esquema MLST','alelo #1','alelo #2','alelo #3','alelo #4','alelo #5','alelo #6','alelo #7','MLSTs posibles','Alelos posibles']
    for col in columnas_finales:
        if col not in resultado_final.columns: resultado_final[col] = ""

    resultado_final = resultado_final[columnas_finales]
    resultado_final.to_excel(taxonomy_out, index=False)
    resultado_final.to_csv(taxonomy_csv, index=False)

    ultimo_df = resultado_final[['Sample','Género mayoritario','Especie mayoritaria','Subespecie','Esquema MLST','MLST','alelo #1','alelo #2','alelo #3','alelo #4','alelo #5','alelo #6','alelo #7','MLSTs posibles','Alelos posibles']].copy()

    if not ultimo_df.empty:
        provis1 = ultimo_df['Género mayoritario'].astype(str).str.split(" ", n=1, expand=True)
        ultimo_df["Género"]     = provis1[0] if 0 in provis1.columns else ""
        ultimo_df["Género (%)"] = provis1[1].str.strip().str.rstrip(')').str.lstrip('(') if 1 in provis1.columns else ""
        provis2 = ultimo_df['Especie mayoritaria'].astype(str).str.rsplit(" ", n=1, expand=True)
        ultimo_df["Especie"]    = provis2[0] if 0 in provis2.columns else ""
        ultimo_df["Especie(%)"] = provis2[1].str.strip().str.rstrip(')').str.lstrip('(') if 1 in provis2.columns else ""
    else:
        ultimo_df["Género"] = ""; ultimo_df["Género (%)"] = ""; ultimo_df["Especie"] = ""; ultimo_df["Especie(%)"] = ""

    ultimo_df.drop(columns=['Género mayoritario', 'Especie mayoritaria'], inplace=True, errors='ignore')
    ultimo_df["Contaminante"] = ""
    ultimo_df = ultimo_df[['Sample',"Género (%)", "Género","Especie(%)", "Especie", "Contaminante",'Subespecie','Esquema MLST','MLST','alelo #1','alelo #2','alelo #3','alelo #4','alelo #5','alelo #6','alelo #7','MLSTs posibles','Alelos posibles']]

    ultimo_df.to_excel(os.path.join(out_folder, "kraken_mlst.xlsx"), index=False)

    print("\n\033[92m[INFO]\033[0m ¡Proceso terminado de forma segura!")

if __name__ == '__main__':
    main()