#!/home/usuario/miniconda3/envs/nano/bin/python

#############################################################
# Jorge R Grande - HUMV - Santander
# parser.py takes in a folder containing the output of minion.sh
# and outputs a short report summarizing stats from the MinION run
# including, QC of reads and assemblies, resistance genes, MLST and plasmids. 



import os
import pandas as pd
import numpy as np
import argparse


def get_arguments():

    """
    Esta función permite introducir los parámetros al llamar al programa en la terminal
    """

    parser = argparse.ArgumentParser(prog = 'parser.py', description = 'Parser.py is part of the ION analysis pipeline. It takes in the resulting csv files from ion.sh and outputs the final report in EXCEL format')

    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input_folder', dest="input_folder", required=True, help="Required. Input folder to be summarized", type=os.path.abspath)

    input_group.add_argument('-db', '--pubMLST_database', dest="pubMLST_database", required=False, default="/home/usuario/miniconda3/envs/mlst/db/pubmlst/", help="Folder containing the MLST schemas", type=os.path.abspath)

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

DB                  = args.pubMLST_database

report              = input_folder + "/mlst.csv"           # Creo que no hace falta ninguna de estas "/"

report_2            = out_folder   + "/mlst_modif.csv"

species_report      = input_folder + "/04_taxonomies/kraken2/species.csv"

genus_report        = input_folder + "/04_taxonomies/kraken2/genus.csv"

gambit_report       = input_folder + "/04_taxonomies/gambit.csv"

kleb_report         = input_folder + "/klebsiellas/enterobacterales__species_output.txt"

ec                  = input_folder + "/04_taxonomies/ectyper/output.tsv"

abricate            = input_folder + "/AbR_report.csv"

abricate_out        = out_folder   + "/AbR_modif.xlsx"

copla_in            = input_folder + "/copla.txt"

copla_out           = out_folder   + "/copla_modif.csv"

taxonomy_out        = out_folder   + "/taxonomy.xlsx"

taxonomy_csv        = out_folder   + "/taxonomy.csv"

########


col_species = ['Sample','Perc. reads', 'Especie']
col_generos = ['Sample','Perc. reads', 'Género']

#kleborate= pd.read_csv(kleb_report, sep='\t', usecols=['strain','species','K_locus','K_locus_confidence','O_locus','O_locus_confidence','Bla_acquired','Bla_ESBL_acquired','Bla_Carb_acquired', 'virulence_score','resistance_score','num_resistance_genes'])

species  = pd.read_csv(species_report, sep= '\t', names=col_species, header=None, decimal = '.')

genus    = pd.read_csv(genus_report, sep= '\t', names=col_generos, header=None, decimal = '.')

gambit   = pd.read_csv(gambit_report)


#################
# Working: MLST #
#################

bad_symb  = ['?', '~']

cabecera = ['Sample','Esquema MLST','MLST','alelo #1','alelo #2','alelo #3','alelo #4','alelo #5','alelo #6','alelo #7','alelo #8','MLSTs posibles','Alelos posibles']

with open(report_2, 'w') as f:
    f.write('\t'.join(map(str, cabecera)) +'\n')
    f.close()

with open(report, 'r') as file:
    for line in file:                                                       # iterar sobre cada muestra:
        column = line.strip('\n').strip().split('\t')                       # eliminar saltos de línea y espacios 


        organismo = column[1]                                               # el mlst.csv indica en la 2ºcol.qué organismo cree que está analizando
        muestra   = column[0].split('/')[1].split('.')[0]
        column[0] = muestra                                                 # parsear el nombre de la muestra
        db_org    = DB + '/' + organismo + '/' + organismo + '.txt'         # coge la DB del organismo

        if (column[1] == '-'):                                              # Si no hay un modelo de MLST escribe la línea de muestra y el error
            column[1] = 'No tiene esquema asociado'
            with open(report_2, 'a') as f:
                f.write('\t'.join(map(str, column)) +'\n')            

        elif (column[1] != '-') and (column[2] != '-'):                     # Si tiene ya MLST escríbela tal cual
            with open(report_2, 'a') as f:
                f.write('\t'.join(map(str, column)) +'\n')          

        elif (column[1] != '-') and (column[2] == '-'):                     # Pero si los ST son indefinidos pero de especie conocida
            genes   = [col.split('(')[0] for col in column[3:]]             # saca la lista provisional de genes de cada muestra, puede incluir genes mal tipados
            numeros = [col.split('(')[1].strip(')') for col in column[3:]]  # saca la lista de alelos de cada muestra

            diccion = dict(zip(genes, numeros))                             # genero un diccionario gen:alelo con ellos

            diccion = {k: v for k, v in diccion.items() if v.isdigit()}     # toma solo los alelos seguros (quita los (?) y (~))
           
            mlst    = pd.DataFrame([diccion])                               # genero un dataframe con el diccionario
            mlst    = mlst.astype('int64')                                  # convierte los alelos a nº entero

            df       = pd.read_csv(db_org, sep = '\t')                      # abro como dataframe la DB de pubMLST

            genes    = list(diccion.keys())                                 # saca una lista de los genes seguros de la muestra
            
            if len(genes) == 0:
                print('No hay ningún alelo bien secuenciado')
                posibles = pd.DataFrame([{'ST': []}])
            
            else:
                posibles = pd.merge(mlst, df, how='left', on=genes)         # mergea la DB de pubMLST con el dataframe de alelos seguros de la muestra

            if (posibles.shape[0] == 1 and df['ST'].isna):
                print("aquí está el problema: no hay ST asociado")

                posibles_alelos = posibles.iloc[:,posibles.columns.get_indexer(['ST'])]

                
            else:
                posibles_alelos = posibles.iloc[:,posibles.columns.get_indexer(['ST'])+1] # aquí me daba el error "positional indexers are out of bounds"

            posibles_alelos = posibles_alelos.to_dict('list')

            for k, v in posibles_alelos.items():
                #v = ['' if pd.isna(x) else x for x in v]
                v = ['Ningún alelo de este gen coincide con esta combinación' if pd.isna(x) else x for x in v]
                posibles_alelos[k] = v[:9]                    

            posibles.fillna('Posible nuevo ST', inplace=True)
            posibles_ST     = posibles.ST.to_list()                         # saca una lista de los posibles ST

            if (len(column) <= 9):
                column.append(str(''))     

            posibles_ST = posibles_ST[:9]
            print(posibles_alelos)

            column.append(posibles_ST)                                      # Ponla al final de la línea del mlst.csv
            column.append(posibles_alelos)                                  # Y finalmente los alelos

            with open(report_2, 'a') as f:

                f.write('\t'.join(map(str, column)) +'\n')
                f.close()

#################
# Working: AbR  #
#################


df_abricate = pd.read_csv(abricate, sep='\t')


for col in df_abricate.columns[2:]:

    
    df_abricate[col] = df_abricate[col].apply(lambda x: col + ' (' + str(x) +')' if x != '.' else '')


diccionario={}

for i in range(len(df_abricate)):

    value = (df_abricate.loc[i, :].values.tolist())

    muestra = value[0]
    valores = [valor for valor in value if valor != ''][2:]
    diccionario[muestra]=valores

df_abricate["Genes resistencia"] = df_abricate["#FILE"].map(diccionario)
df_abricate["Genes resistencia"] = df_abricate["Genes resistencia"].apply(lambda x: ', '.join([i for i in x]))
df_abricate["#FILE"] = df_abricate["#FILE"].str.split('/').str[1]

cols = df_abricate.columns.tolist()
cols = cols[:2] + cols[-1:] + cols[2:-1]
df_abricate  = df_abricate[cols]

df_abricate.to_excel(abricate_out, index=False)


#################
# Working: copla#
#################


d={'Sample': [], 'Contig': [], 'Size': [], 'MOB': [], 'MPF': [], 'Rep': [], 'AbR': []}


with open(copla_in, 'r') as f:
    while True:
        line=f.readline()
        line=line.lstrip('  ')
        if line.startswith('Sample'):
                print(line)
                if line.split(':')[1].strip() != '':
                     d['Sample'].append(line.split(':')[1].strip())
                else:
                     d['Sample'].append('-')
        if line.startswith('Contig'):
                print(line)
                if line.split(':')[1].strip() != '':
                     d['Contig'].append(line.split(':')[1].strip())
                else:
                     d['Sample'].append('-')              
        if line.startswith('Size'):
                print(line)
                if line.split(':')[1].strip() != '':
                     d['Size'].append(line.split(':')[1].strip())
                else:
                     d['Sample'].append('-')
        if line.startswith('MOB'):
                print(line)
                if line.split(':')[1].strip() != '':
                     d['MOB'].append(line.split(':')[1].strip())
                else:
                     d['Sample'].append('-')
        if line.startswith('MPF'):
                print(line)
                if line.split(':')[1].strip() != '':
                     d['MPF'].append(line.split(':')[1].strip())
                else:
                     d['Sample'].append('-')
        if line.startswith('Repl'):
                print(line)
                if line.split(':')[1].strip() != '':
                     d['Rep'].append(line.split(':')[1].strip())
                else:
                     d['Sample'].append('-')
        if line.startswith('AMR'):
                print(line)
                if line.split(':')[1].strip() != '':
                     d['AbR'].append(line.split(':')[1].strip())
                else:
                     d['Sample'].append('-')
        if not line:
            break

df = pd.DataFrame(d).set_index('Sample')

df.to_csv(copla_out)

##################################
## Taxonomia - MLST - Kleborate ##
##################################

mlst     = pd.read_csv(report_2, sep='\t')
mlst['Sample'] = mlst['Sample'].astype(str)
print("MLST ", mlst)
print(mlst.dtypes)

col_species = ['Sample','Per reads', 'Especie']
species  = pd.read_csv(species_report, sep= '\t', names=col_species, header=None, decimal = '.')

col_generos = ['Sample','Per reads', 'Género']
genus    = pd.read_csv(genus_report, sep= '\t', names=col_generos, header=None, decimal = '.')

species['Sample'] = species['Sample'].astype(str)
genus['Sample']   = genus['Sample'].astype(str)

gambit   = pd.read_csv(gambit_report)

kleborate= pd.read_csv(kleb_report, sep='\t', usecols=['strain','enterobacterales__species__species','klebsiella_pneumo_complex__kaptive__K_locus',
							'klebsiella_pneumo_complex__kaptive__K_locus_confidence',
							'klebsiella_pneumo_complex__kaptive__O_locus',
							'klebsiella_pneumo_complex__kaptive__O_locus_confidence',
							'klebsiella_pneumo_complex__amr__Bla_acquired',
							'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired',
							'klebsiella_pneumo_complex__amr__Bla_Carb_acquired',
							'klebsiella_pneumo_complex__virulence_score__virulence_score',
							'klebsiella_pneumo_complex__resistance_score__resistance_score',
							'klebsiella_pneumo_complex__resistance_gene_count__num_resistance_genes'])
kleborate['strain'] = kleborate['strain'].astype(str)
kleborate.fillna('', inplace=True)
gambit = gambit[['query','closest.description']]
gambit = gambit.rename(columns={'query':'Sample','closest.description':'Subespecie'})
gambit['Subespecie'] = gambit['Subespecie'].str.replace(r"\(.*?\)", "", regex=True)
gambit['Subespecie'] = gambit['Subespecie'].str.replace(r"\[.*?\]", "", regex=True)
gambit['Sample'] = gambit['Sample'].astype(str)
gambit = gambit.applymap(lambda x: x.strip() if isinstance(x, str) else x)


species['Especie mayoritaria'] = species['Especie'] + ' (' + species["Per reads"].round(0).astype(int).astype(str) +  '%)'
genus  ['Género mayoritario']  = genus['Género'] + ' (' + genus['Per reads'].round(0).astype(int).astype(str) +  '%)'
genus.drop(['Per reads','Género'], axis = 1, inplace=True)

species['Posibles contaminantes'] = species.groupby(['Sample'])['Especie mayoritaria'].transform(lambda x: '-' if len(x) == 1 else ','.join(x.iloc[1:]))
species = species.drop_duplicates(subset=['Sample','Posibles contaminantes'], keep = 'first')
species.drop(['Per reads','Especie'], axis = 1, inplace=True)

kraken  = pd.merge(species, genus, how = "left", on='Sample')

kraken = kraken[['Sample','Género mayoritario','Especie mayoritaria','Posibles contaminantes']]
kraken['Sample'] = kraken['Sample'].astype(str)
kraken = kraken.applymap(lambda x: x.strip() if isinstance(x, str) else x)

subesp = pd.merge(kraken, gambit, how="left", on='Sample')

intermedio = pd.merge(subesp, mlst, how='left', on='Sample')

kleborate['klebsiella_pneumo_complex__amr__Bla_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_acquired'].str.replace(r'\.v1\^', '', regex=True)
kleborate['klebsiella_pneumo_complex__amr__Bla_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_acquired'].str.replace(r'\^', '', regex=True)
kleborate['klebsiella_pneumo_complex__amr__Bla_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_acquired'].str.split(';')
kleborate['klebsiella_pneumo_complex__amr__Bla_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_acquired'].apply(set).apply(list).apply(lambda x: ', '.join(map(str, x)))

kleborate['klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'].str.replace(r'\.v1\^', '', regex=True)
kleborate['klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'].str.replace(r'\^', '', regex=True)
kleborate['klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'].str.split(';')
kleborate['klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_ESBL_acquired'].apply(set).apply(list).apply(lambda x: ', '.join(map(str, x)))

kleborate['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'].str.replace(r'\.v1\^', '', regex=True)
kleborate['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'].str.replace(r'\^', '', regex=True)
kleborate['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'].str.split(';')
kleborate['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'] = kleborate['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'].apply(set).apply(list).apply(lambda x: ', '.join(map(str, x)))

kleborate['klebsiella_pneumo_complex__kaptive__K_locus']    = kleborate['klebsiella_pneumo_complex__kaptive__K_locus'].str.replace(r'unknown \([A-Z]*[0-9]*\-*[A-Z]*[0-9]*\)', '-', regex=True)  
kleborate['klebsiella_pneumo_complex__kaptive__O_locus']    = kleborate['klebsiella_pneumo_complex__kaptive__O_locus'].str.replace(r'unknown \([A-Z]*[0-9]*\/*[A-Z]*[0-9]*[av]*[0-9]*\)', '-', regex=True)
kleborate['K/O locus']  = kleborate['klebsiella_pneumo_complex__kaptive__K_locus']            + '/' + kleborate['klebsiella_pneumo_complex__kaptive__O_locus']

kleborate.drop(['klebsiella_pneumo_complex__kaptive__K_locus_confidence','klebsiella_pneumo_complex__kaptive__O_locus_confidence','klebsiella_pneumo_complex__kaptive__O_locus','enterobacterales__species__species'], axis = 1, inplace=True)
kleborate = kleborate.rename(columns={'strain':'Sample',
					'klebsiella_pneumo_complex__amr__Bla_acquired':'Otras',
					'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired':'BLEE adquirida',
					'klebsiella_pneumo_complex__amr__Bla_Carb_acquired':'Carba adquirida',
					'klebsiella_pneumo_complex__virulence_score__virulence_score':'VIRscore',
					'klebsiella_pneumo_complex__resistance_score__resistance_score':'AMRscore',
					'klebsiella_pneumo_complex__resistance_gene_count__num_resistance_genes':'Nº genes AMR'})

kleborate = kleborate[['Sample', 'K/O locus','Carba adquirida','BLEE adquirida','Otras', 'Nº genes AMR','AMRscore','VIRscore']]

resultado = pd.merge(intermedio, kleborate, how="left", on='Sample')

# 

ectyper  = pd.read_csv(ec, sep='\t', usecols=['Name','Serotype'])
ectyper  = ectyper.rename(columns={'Name':'Sample'})
ectyper['Sample'] = ectyper['Sample'].astype(str)


resultado_final  = pd.merge(resultado, ectyper, how="left", on='Sample')
resultado_final  = resultado_final[['Sample','Género mayoritario','Especie mayoritaria','Subespecie','MLST','Serotype','K/O locus','Posibles contaminantes','Carba adquirida','BLEE adquirida','Otras','Nº genes AMR','AMRscore','VIRscore','Esquema MLST','alelo #1','alelo #2','alelo #3','alelo #4','alelo #5','alelo #6','alelo #7','MLSTs posibles','Alelos posibles']]
# resultado_final.drop(columns=['confidence'], axis = 1, inplace=True)

resultado_final.to_excel(taxonomy_out, index=False)
resultado_final.to_csv  (taxonomy_csv, index=False)
