# %%
import pandas as pd
from glob import glob
import dateutil.parser
import re
import os
import argparse
import numpy as np

# %%
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script para procesar datos de secuenciación y ensamblaje")
    
    # Argumento obligatorio para 'base_run'
    parser.add_argument("--input_path", type=str, help="Path del run de secuenciación")

    # Argumento opcional para 'output_run'
    parser.add_argument("--output_file", type=str, default="datos_seq_nuevos.tsv", 
                        help="Nombre del archivo de salida (opcional, por defecto 'output_results.tsv')")

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
        d_sum['fecha'] = d_sum['start'].strftime('%Y-%m-%d')
        hours, minutes, seconds = str(duration).split(sep=':')
        d_sum['duracion'] = hours + 'h ' + minutes + 'min.'

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
            # Parseo de los links
            if line.startswith('L'):
                edge1 = line.split()[1]
                edge2 = line.split()[3]
                if edge1 not in d_links:
                    d_links[edge1] = []
                d_links[edge1].append(edge2)
                d_links[edge1] = list(set(d_links[edge1]))
            # Parseo de los paths
            if line.startswith('P'):
                contig = line.split()[1]
                edges = line.split()[2].split(',')
                edges = [re.sub(r"[+-]", "", edge) for edge in edges]
                d_paths[contig] = set(edges)
    return d_links, d_paths

# %%
def parse_info(file_path, d_links, d_paths):
    df = pd.read_table(file_path)
    df.sort_values('length', ascending=False, inplace=True)    # Por defecto está ordenado así, pero por si acaso
    d_info = {}

    # contig mayor cerrado?
    d_info['contig1_closed'] = df['circ.'][0] == 'Y'

    # ratio de tamaño (2 contigs mayores / total)
    contig1 = df['#seq_name'][0]
    if len(df) != 1:
        total_length = sum(df['length'])
        l_contig1 = df['length'][0]
        l_contig2 = df['length'][1]
        d_info['ratio_total'] = (l_contig1 + l_contig2) / total_length

        # ratio entre los dos contigs mayores
        d_info['ratio_greatest_contigs'] = l_contig2 / l_contig1

        # Comprobar si los dos contigs mayores están unidos por un nodo común
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

    # número de plásmidos
    if d_info['are_linked']:
        chromosomes = {contig1, contig2}
    else:
        chromosomes = {contig1}
    plasmids = [cont for cont in df['#seq_name'] if cont not in chromosomes]
    d_info['n_plas'] = len(plasmids)

    # número de plásmidos cerrados
    closed_plas = df[df['#seq_name'].isin(plasmids)]['circ.']
    closed_plas.reset_index()
    d_info['n_closed_plas'] = sum(closed_plas == 'Y')

    # return n_contigs, n_closed, contig1_closed, ratio_total, ratio_greatest_contigs, are_linked, contig1, contig2, n_plas, closed_plas
    return d_info

# %%
# Cálculo de puntuación
def calculate_score(d_info):

    score = 0
    d_score = {}

    # Puntuación con respecto al cromosoma
    # ratio de los dos mayores contigs/total
    if d_info['ratio_total'] > 0.8: 
        score += 1.5
        d_score['ratio'] = 1.5
    else:
        d_score['ratio'] = 0

    # Si el segundo contig más largo es como máximo un 10% del contig más largo
    if d_info['ratio_greatest_contigs'] < 0.1:
        score += 1.25
        d_score['ratio_1_2'] = 1.25
        # Si el cromosoma está cerrado
        if d_info['contig1_closed']:
            score += 0.75
            d_score['contig1_closed'] = 0.75
        else:
            d_score['contig1_closed'] = 0
    # Si no, valorar si los contigs están unidos (secuencia de inserción, 0.5) o no
    else:
        d_score['ratio_1_2'] = 0
        if d_info['are_linked'] == True:
            score += 0.75
            d_score['greatest_linked'] = 0.75
        else:
            d_score['greatest_linked'] = 0
    
    # Puntuación con respecto a plásmidos
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
    for sample in lista_cepas['ID único']:
        assembly_path = f'{base_run}/03_assemblies/{sample}/'
        info = os.path.join(assembly_path, 'assembly_info.txt')
        gfa = os.path.join(assembly_path, 'assembly_graph.gfa')

        # Si el assembly se ha ensamblado calcula la puntuación, si no devuelve puntuación 0
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

def main():
    # Input path con todos los archivos
    args = parse_arguments()
    base_run = args.input_path
    output_run = args.output_file

    # %%
    # Input files
    # Resúmenes finales de la secuenciación (MinION)
    summary = glob(os.path.join(base_run, "final_summary*txt"))[0]  # glob saca una lista de ficheros que se parecen a lo que le pides, no te saca el string con el path directo!
    report = glob(os.path.join(base_run, "report_*.json"))[0]

    # Estadísticas control de calidad
    qc_r = os.path.join(base_run, "QC_reads.csv")
    qc_a = os.path.join(base_run, "QC_assembly.csv")

    # Histórico de datos de secuenciación
    tabla = os.path.join(base_run, "datos_seq.tsv")

    # Histórico de datos de análisis
    anali = os.path.join(base_run, "datos_analisis.tsv")

    # Información básica del run actual 
    cepas = os.path.join(base_run, "lista_seq.tsv")

    # Taxonomy
    taxon = os.path.join(base_run, "taxonomy.csv")
    
    # EGMs
    # Plásmidos
    plasmids = os.path.join(base_run, "copla_modif.tsv")

    # ICEs
    ices = os.path.join(base_run, "ICE_summary.csv")

    # Fagos
    fagos = os.path.join(base_run, "phage_summary.csv")

    # Integrones
    integrones = os.path.join(base_run, "integron_summary.csv")

    # Comprobamos que todos los archuvos necesarios existen
    required_files = [summary, report, qc_r, qc_a, tabla, cepas, taxon]
    missing_files = [f for f in required_files if not os.path.isfile(f)]
    if missing_files:
        raise FileNotFoundError(f"Los archivos siguientes no existen: {missing_files}")


    # %%
    # Cargamos los inputs
    # Tabla con las muestras secuenciadas y pendientes de secuenciación
    datos_seq = pd.read_csv(tabla, sep='\t')
    # Borramos todas las columnas vacías (Unnamed: )
    datos_seq = datos_seq.loc[:, ~datos_seq.columns.str.contains('^Unnamed: ')]
    # Definir formato de columnas
    datos_seq['Nº Cultivo'] = datos_seq['Nº Cultivo'].astype(str)
    datos_seq["Barcode"]    = datos_seq["Barcode"].astype(str)
    datos_seq['Barcode'] = datos_seq['Barcode'].replace('nan', np.nan)

    # Tabla con información del run actual
    lista_cepas = pd.read_csv(cepas, sep='\t', usecols=["Nº Cultivo", "Cepario", "ID único", "Barcode", "[DNA]"], dtype={'Barcode': 'string'})
    lista_cepas['Barcode'] = lista_cepas['Barcode'].str.replace(r'barcode', '', regex=True)

    # Añadir muestras nuevas
    nuevas_filas = lista_cepas[~lista_cepas['ID único'].isin(datos_seq['ID único'])]
    datos_seq = pd.concat([datos_seq, nuevas_filas], ignore_index=True)

    # Diccionarios con información técnica de la secuenciación
    d_sum = parse_minion_sum(summary)
    d_report = parse_minion_report(report)

    # Tablas con calidad de la secuenciación y los assemblies
    QC_reads = pd.read_csv(qc_r, sep='\t', decimal='.', thousands=',')
    QC_assembly = pd.read_csv(qc_a, sep='\t')

    # Diccionario con información de los ensamblajes
    d_quality = get_assembly_score(lista_cepas, base_run)

    # %%
    # Output files
    # Output datos secuenciación
    output_run = os.path.join(base_run, output_run)
    analisis_run = os.path.join(base_run,"datos_analisis_nuevos.tsv")


    # Inicializamos la tabla output con las muestras de este run y la información técnica
    columnas = ["Nº Cultivo", "Cepario", "ID único", "Fecha seq", "Fecha seq (rep)", "Fecha seq (rep2)", 
                "Kit Extracción", "Kit Barcoding", "Barcode", "BC (rep)", "BC (rep2)", "Posición", 
                "Tipo FlowCell", "FlowCell", "Poros inicio", "Poros final", "Horas seq", "Cepas/run", 
                "Cepas a repetir/run", "Rendim (Mbp)", "Repetir", "Tª trabajo", "Voltaje ", "Rend/h (reads)", 
                "Rend/h (MbP)", "N50 (kbp)"]
    
    # Creamos 'result' copiando 'lista_cepas' y añadimos las columnas que falten como nulas
    result = lista_cepas.copy()
    for col in columnas:
        if col not in result.columns:
            result[col] = pd.NA

    # %%
    # Actualizamos la tabla con la información técnica
    result["Fecha seq"] = d_sum['fecha']
    result["Kit Barcoding"] = d_sum['barcoding_kit']
    result["Kit Extracción"] = "DNeasy Blood & Tissue"
    result["Posición"] = d_sum['instrument']
    result["Tipo FlowCell"] = d_sum['flow_cell_type']
    result["FlowCell"] = d_sum['flow_cell']
    result["Horas seq"] = d_sum['duracion']
    result["Poros inicio"] = d_report['poros_ini']
    result["Poros final"] = d_report['poros_fin']

    # %%
    # Poblar con datos de QC_reads.csv
    QC_reads = QC_reads.rename(columns={"Sample":"ID único",
                            "Median length" : "Lmediana (pre)", 
                            "Median quality" : "Qmediana (pre)", 
                            "Total reads" : "Nreads (pre)", 
                            "Total bases" : "Nbases (pre)",
                            "Median length.1" : "Lmediana (post)", 
                            "Median quality.1" :"Qmediana (post)", 
                            "Total reads.1" : "Nreads (post)",
                            "Total bases.1" : "Nbases (post)"})
                            
    QC_reads = QC_reads.drop(columns=["MaxQ", "Longest read", "Sample.1", "MaxQ.1", "Longest read.1"])
    QC_reads[["Nbases (post)"]].apply(pd.to_numeric)
    result2 = pd.merge(result, QC_reads, on="ID único", how='outer')

    # Poblar con datos de QC_assembly.csv
    QC_assembly = QC_assembly.rename(columns={"Sample":"ID único"})
    QC_assembly["ratio"] = QC_assembly["Largest contig"]/QC_assembly["Total length"]
    QC_assembly = QC_assembly.drop(columns=["GC (%)",	"# predicted genes (>= 300 bp)"], errors='ignore')

    result3 = pd.merge(result2, QC_assembly, on="ID único", how='outer')

    result3[["Total length", "Nbases (post)"]].apply(pd.to_numeric)

    result3["Profundidad"] = result3["Nbases (post)"].div(result3["Total length"])
    result3["Profundidad"] = result3["Profundidad"].round(0).astype('Int64')

    result3["% Bases Filtrado"] = result3["Nbases (post)"].div(result3["Nbases (pre)"])
    
    # Añadir score del assembly
    result3['Calidad (ass)'] = result3['ID único'].map(d_quality)



    orden_final = ["Nº Cultivo", "Cepario","ID único", "Barcode", "BC (rep)", "BC (rep2)", "Fecha seq", "Fecha seq (rep)", "Fecha seq (rep2)", "[DNA]",
            "Profundidad", 'Calidad (ass)', "Kit Extracción", "Kit Barcoding", "Posición", "Tipo FlowCell", 
            "FlowCell", "Poros inicio", "Poros final", "Horas seq", "Cepas/run", "Cepas a repetir/run", 
            "Lmediana (pre)", "Qmediana (pre)", "Nreads (pre)", "Nbases (pre)",          
            "Lmediana (post)", "Qmediana (post)", "Nreads (post)", "Nbases (post)", "% Bases Filtrado"]
    result3 = result3[orden_final]

    # %%
    Ncepas_inicial = lista_cepas.shape[0]
    Ncepas_bien = (result3["Profundidad"] > 30.0).sum()
    Ncepas_repetir = Ncepas_inicial - Ncepas_bien
    
    result3["Cepas/run"]           = Ncepas_inicial
    result3["Cepas a repetir/run"] = Ncepas_repetir

    result3['Nº Cultivo'] = result3['Nº Cultivo'].astype(str)
    result3["Barcode"]    = result3["Barcode"].astype(str)

    # Realiza el merge para unir las tablas basado en 'ID único'
    merged_df = pd.merge(datos_seq, result3, on='ID único', how='left', suffixes=('', '_result3'))
    
    # Verifica y asigna los valores de "Fecha seq (rep2)" y "BC (rep2)"
    merged_df['Fecha seq (rep2)'] = merged_df['Fecha seq (rep2)'].combine_first(
        merged_df.apply(lambda x: x['Fecha seq_result3'] if pd.notna(x['Fecha seq (rep)']) else None, axis=1))
    merged_df['BC (rep2)'] = merged_df['BC (rep2)'].combine_first(
        merged_df.apply(lambda x: x['Barcode_result3'] if pd.notna(x['BC (rep)']) else None, axis=1))

    # Verifica y asigna los valores de "Fecha seq" y "Barcode"
    merged_df['Fecha seq (rep)'] = merged_df['Fecha seq (rep)'].combine_first(
        merged_df.apply(lambda x: x['Fecha seq_result3'] if pd.notna(x['Fecha seq']) else None, axis=1))
    merged_df['BC (rep)'] = merged_df['BC (rep)'].combine_first(
        merged_df.apply(lambda x: x['Barcode_result3'] if pd.notna(x['Barcode']) else None, axis=1))

    # Finalmente, rellena cualquier valor faltante en las columnas originales
    merged_df['Fecha seq'] = merged_df['Fecha seq'].combine_first(merged_df['Fecha seq_result3'])
    merged_df['Barcode'] = merged_df['Barcode'].combine_first(merged_df['Barcode_result3'])






    # %%
    # Rellena las filas vacías de datos_seq con las de result3
    for column in datos_seq.columns:
        if column not in ['ID único', 'Barcode', 'BC (rep)', 'BC (rep2)', 'Fecha seq', 'Fecha seq (rep)', 'Fecha seq (rep2)']:  # Evitar la columna de unión
            merged_df[column] = merged_df[column + '_result3'].combine_first(merged_df[column])

    # Elimina las columnas extra de result3
    merged_df = merged_df[datos_seq.columns]

    # Redondeo las columnas con valores enteros para que no queden con el ".0 tan feo
    int_columns = [
        'Profundidad','Cepas/run','Cepas a repetir/run','Lmediana (pre)',
        'Nreads (pre)','Nbases (pre)','Lmediana (post)','Nreads (post)','Nbases (post)'
    ]
    for col in int_columns:
        if col in merged_df.columns:
            merged_df[col] = merged_df[col].astype(str).str.replace(',', '', regex=False).str.strip()
            merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce').round().astype('Int64')
    # Y redondeo la última columna a solo 2 decimales:
    merged_df["% Bases Filtrado"] = (merged_df["% Bases Filtrado"]*100).round(1)
    merged_df.to_csv(output_run, index=False, sep = "\t")

    # %%
    analisis = pd.read_csv(anali, sep='\t')
    analisis.rename(columns={"Muestra":"Nº Cultivo", "Serotipo": "Serotype"}, inplace=True)

    taxon2 = pd.read_csv(taxon, sep= ',')
    taxon2.rename(columns={"Sample":"ID único"}, inplace=True)
    taxon2.drop_duplicates(subset=['ID único'], keep='first', inplace=True)

    result4 = pd.merge(taxon2, lista_cepas, on="ID único", how='outer')
    result4.drop(columns=["Cepario", "[DNA]"], inplace=True)



    # EGMs
    result4[["Plásmidos", "ICEs", "Profagos", "Integrones"]] = 0

    # Plásmidos
    plasmids = os.path.join(base_run, "copla_modif.tsv")
    df_pl = pd.read_csv(plasmids, sep='\t')
    pl_count = df_pl['Sample'].value_counts()
    result4['Plásmidos'] = result4['ID único'].map(pl_count, na_action='ignore')

    # ICEs
    # ices = os.path.join(base_run, "ICE_summary.csv")
    # df_ices = pd.read_csv(ices, sep=',')
    # ice_count = df_ices['Nombre muestra'].value_counts()
    # result4['ICEs'] = result4['ID único'].map(ice_count, na_action='ignore')
    result4["ICEs"] = '0'

    # Fagos
    fagos = os.path.join(base_run, "phage_summary.csv")
    df_fagos = pd.read_csv(fagos, sep=',')
    fago_count = df_fagos['sample'].value_counts()
    result4['Profagos'] = result4['ID único'].map(fago_count, na_action='ignore')

    # Integrones
    integrones = os.path.join(base_run, "integron_summary.csv")
    df_int = pd.read_csv(integrones, sep=',')
    int_count = df_int['Sample'].value_counts()
    result4['Integrones'] = result4['ID único'].map(int_count, na_action='ignore')

    result4[['Plásmidos', 'ICEs', 'Profagos', 'Integrones']] = result4[['Plásmidos', 'ICEs', 'Profagos', 'Integrones']].fillna(0)

    result4 = result4.merge(merged_df[['ID único', 'Calidad (ass)', 'Profundidad']], on='ID único', how='inner')

    nwo = ["Nº Cultivo", "ID único", "Barcode",  "Profundidad", 'Calidad (ass)',	
       "Subespecie", "MLST", "Serotype", "K/O locus", "Posibles contaminantes",	"Carba adquirida",	
       "BLEE adquirida", "Otras", "Nº genes AMR", "AMRscore", "VIRscore", 
       "Plásmidos",  "Profagos", "Integrones", 	
       "alelo #1",	"alelo #2",	"alelo #3",	"alelo #4",	"alelo #5",	"alelo #6",	"alelo #7",	
       "MLSTs posibles", "Alelos posibles",  "Género mayoritario", "Especie mayoritaria","Esquema MLST"]
    result4=result4[nwo]
    analisis=analisis[nwo]

    analisis['ID único'] = analisis['ID único'].astype(str)
    analisis["Barcode"]    = analisis["Barcode"].astype(str)
    result4['ID único'] = result4['ID único'].astype(str)
    result4["Barcode"]    = result4["Barcode"].astype(str)
    analisis['Barcode'] = analisis['Barcode'].replace('nan', np.nan)




    # Realiza el merge para unir las tablas basado en 'ID único'
    analisis_final = pd.merge(analisis, result4, on='ID único', how='left', suffixes=('', '_result4'))

    for column in analisis.columns:
        if column != 'ID único':  # Evitar la columna de unión
            analisis_final[column] = analisis_final[column + '_result4'].combine_first(analisis_final[column])    # está fallando porque no tienen el mismo orden las dos tablas

    # Elimina las columnas extra de result4
    analisis_final = analisis_final[analisis.columns]

    # Sustituyo barcode13 por solo 13
    analisis_final['Barcode'] = analisis_final['Barcode'].str.replace(r'barcode', '', regex=True)

    # Redondeo las columnas con valores enteros para que no queden con el ".0 tan feo
    int_columns = [
        "Nº Cultivo",'Profundidad',"Nº genes AMR","AMRscore", "VIRscore", "Plásmidos", "Profagos", "Integrones"
    ]
    for col in int_columns:
        if col in analisis_final.columns:
            analisis_final[col] = analisis_final[col].astype(str).str.replace(',', '', regex=False).str.strip()
            analisis_final[col] = pd.to_numeric(analisis_final[col], errors='coerce').round().astype('Int64')

    analisis_final.to_csv(analisis_run, index=False, sep='\t')

# %%
if __name__ == "__main__":
    main()
