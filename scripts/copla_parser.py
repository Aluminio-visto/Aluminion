#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd



def run_parsing(input_folder, out_folder=None):
    """Parses the output of COPLA, that usualli looks like this:
    Sample: barcode61
    Contig: AC507
    Query is a PTU-Q2 plasmid
    Query is part of a sHSBM cluster of size 5
    Other info:
    Size: 8254
    MOB:  MOBP
    MPF:  -
    Repl: IncQ2
    AMR:  GES-6
    Sample: barcode61
    Contig: AA002
    Query is a PTU-L/M plasmid
    Query is part of a sHSBM cluster of size 60
    Other info:
    Size: 61273
    MOB:  MOBP
    MPF:  typeI
    Repl: IncL/M(pMU407)
    AMR:  -  

    and creates a tabular output like this one: 
    Sample  Contig  PTU     Size    MOB     MPF     Rep     AbR
    barcode61       AC507   PTU-Q2  16508   MOBP       -       IncQ2     GES-6
    barcode61       AA002   PTU-L/M 61273   MOBP    typeI   IncL/M(pMU407)  -
    
    """
    if out_folder is None:
        out_folder = input_folder

    copla_in = os.path.join(input_folder, "copla.txt")
    copla_out = os.path.join(out_folder, "copla_modif.csv")
    
    records = []         
    current_record = {}  

    if os.path.exists(copla_in):
        with open(copla_in, 'r') as f:
            for line in f:
                line = line.strip() 
                if not line or line == "Other info:":
                    continue
                    
                if line.startswith('Sample:'):
                    if current_record:
                        records.append(current_record)
                        current_record = {} 
                    current_record['Sample'] = line.split(':')[1].strip()
                elif line.startswith('Contig:'):
                    current_record['Contig'] = line.split(':')[1].strip()
                elif line.startswith('Query is a'):
                    partes = line.split(' ')
                    if len(partes) >= 4:
                        current_record['PTU'] = partes[3].strip()
                elif line.startswith('PTU could not'):
                    current_record['PTU'] = '-'
                elif line.startswith('Size:'):
                    current_record['Size'] = line.split(':')[1].strip()
                elif line.startswith('MOB:'):
                    current_record['MOB'] = line.split(':')[1].strip()
                elif line.startswith('MPF:'):
                    current_record['MPF'] = line.split(':')[1].strip()
                elif line.startswith('Repl:'):
                    current_record['Rep'] = line.split(':')[1].strip()
                elif line.startswith('AMR:'):
                    current_record['AbR'] = line.split(':')[1].strip()

        if current_record:
            records.append(current_record)

        df_copla = pd.DataFrame(records)
        df_copla.fillna('-', inplace=True)
        # Save to CSV so aluminion_reporter.py can find it
        df_copla.to_csv(copla_out, index=False)
        print(f" -> Parsed plasmids saved to {copla_out}")
    else:
        print(f"\033[93m[WARNING]\033[0m COPLA file not found in: {copla_in}")

if __name__ == "__main__":
    import sys
    run_parsing(sys.argv[1])