#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import base64
import os
import sys
import ast
import re

def get_base64_image(image_path):
    """Reads a PNG image and converts it to a Base64 string."""
    if os.path.exists(image_path):
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    return None

def load_and_standardize(file_path, sep=',', key_col=None):
    """Loads a file, standardizes the key column, and returns a DataFrame."""
    if not os.path.exists(file_path):
        print(f"Warning: File not found: {file_path}")
        return pd.DataFrame()

    try:
        if file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path)
        else:
            df = pd.read_csv(file_path, sep=sep)

        # Rename requested key columns (accepts a string or list of candidates)
        rename_map = {}
        candidates = [key_col] if isinstance(key_col, str) else (key_col or [])
        for col in candidates:
            if col and col in df.columns:
                rename_map[col] = 'Sample'
                break
        if 'Cultivo' in df.columns:
            rename_map['Cultivo'] = 'ID'
        if 'Barcode' in df.columns:
            rename_map['Barcode'] = 'BC'

        df.rename(columns=rename_map, inplace=True)

        if 'Sample' in df.columns:
            df['Sample'] = df['Sample'].astype(str).str.strip()

        return df
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return pd.DataFrame()

# ==========================================
# Unit formatting functions (Genomics)
# ==========================================
def fmt_mbp_0(x):
    try: return f"{float(x)/1e6:.0f} Mbp"
    except: return x

def fmt_mbp_2(x):
    try: return f"{float(x)/1e6:.2f} Mbp"
    except: return x

def fmt_kbp_dynamic(x):
    """Uses 1 decimal if below 10 Kbp, and 0 decimals if 10 Kbp or above."""
    try: 
        val = float(x)
        if val < 10000:
            return f"{val/1e3:.1f} Kbp"
        else:
            return f"{val/1e3:.0f} Kbp"
    except: return x

def main():
    work_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    os.chdir(work_dir)
    print(f"Generando informe en: {os.getcwd()}")

    # 1. Load tables
    df_lista    = load_and_standardize("list_seq.tsv", sep='\t', key_col='ID')
    df_qcreads  = load_and_standardize("QC_reads.csv", sep='\t', key_col='Sample')
    df_qcass    = load_and_standardize("QC_assembly.csv", sep='\t', key_col='Samples')
    df_tax      = load_and_standardize("taxonomy.xlsx", key_col='Sample')
    df_abr      = load_and_standardize("AbR_modif.xlsx", key_col='#FILE')
    df_copla    = load_and_standardize("copla_modif.csv", sep=',', key_col='Sample')
    df_integron = load_and_standardize("integron_summary.csv", sep=',', key_col='Sample')
    df_phage    = load_and_standardize("phage_summary.csv", sep=',', key_col='Sample')
    df_kleb     = load_and_standardize("kleborate.tsv", sep='\t', key_col='strain')

    # ==========================================
    # INTEGRON PROCESSING
    # ==========================================
    if not df_integron.empty:
        cassette_cols = [col for col in df_integron.columns if str(col).startswith('Cassette')]
        
        def fix_cassettes(row):
            flattened = []
            for col in cassette_cols:
                val = row[col]
                # Skip empty or null values
                if pd.isna(val) or str(val).strip() == "":
                    continue

                # 1. Aggressively clean hypothetical proteins, brackets and quotes from the string
                clean_val = str(val).replace('[', '').replace(']', '').replace("'", "").replace('"', '').replace(';hypothetical protein', '').replace(';Multidrug transporter EmrE', '')

                # 2. Split by comma (now that the string is clean)
                genes_in_cassette = [g.strip() for g in clean_val.split(',')]

                # 3. Filter out any remaining empty elements
                genes_in_cassette = [g for g in genes_in_cassette if g]
                flattened.extend(genes_in_cassette)
            
            # Asignamos los genes aplanados a las columnas correspondientes
            for i, col in enumerate(cassette_cols):
                row[col] = flattened[i] if i < len(flattened) else ""
            return row

        df_integron = df_integron.apply(fix_cassettes, axis=1)
        
        # Drop Cassette columns above 10 if empty
        cassettes_to_drop = [col for col in cassette_cols if int(col.replace('Cassette', '').strip()) > 10]
        df_integron.drop(columns=cassettes_to_drop, errors='ignore', inplace=True)

    # ==========================================
    # QC TABLE CONSTRUCTION
    # ==========================================
    df_qc = df_lista.copy()
    
    if not df_qcreads.empty:
        rename_dict = {
            'Median length': 'Median L_pre', 'Median quality': 'Median Q_pre',
            'Total reads': 'Total reads_pre', 'Total bases': 'Total bases_pre',
            'MaxQ': 'Max Q_pre', 'Longest read': 'Longest read_pre',
            'Median length.1': 'Median L_post', 'Median quality.1': 'Median Q_post',
            'Total reads.1': 'Total reads_post', 'Total bases.1': 'Total bases_post',
            'MaxQ.1': 'Max Q_post', 'Longest read.1': 'Longest read_post'
        }
        df_qcreads.rename(columns=rename_dict, inplace=True)
        df_qc = pd.merge(df_qc, df_qcreads, on='Sample', how='left')

    if not df_qcass.empty:
        df_qc = pd.merge(df_qc, df_qcass, on='Sample', how='left')
    
    df_qc.drop(columns=['Strain', 'DNA_conc', 'is_repeated', 'Sample.1', '# predicted genes (>= 300 bp)'], errors='ignore', inplace=True)

    int_columns = [
        'Total reads_pre', 'Total bases_pre', 'Longest read_pre', 'Median L_pre',
        'Total reads_post', 'Total bases_post', 'Longest read_post', 'Median L_post',
        '# contigs', 'Largest contig', 'Total length'
    ]
    for col in int_columns:
        if col in df_qc.columns:
            df_qc[col] = df_qc[col].astype(str).str.replace(',', '', regex=False).str.strip()
            df_qc[col] = pd.to_numeric(df_qc[col], errors='coerce').round().astype('Int64')

    if 'Total bases_post' in df_qc.columns and 'Total length' in df_qc.columns:
        tb = df_qc['Total bases_post'].astype(float)
        tl = df_qc['Total length'].astype(float)
        df_qc['Depth'] = (tb / tl).round(0).astype('Int64')

    for c in ['Total bases_pre', 'Total bases_post']:
        if c in df_qc.columns: df_qc[c] = df_qc[c].apply(fmt_mbp_0)
        
    for c in ['Total length', 'Largest contig']:
        if c in df_qc.columns: df_qc[c] = df_qc[c].apply(fmt_mbp_2)
        
    for c in ['Median L_pre', 'Median L_post', 'Longest read_pre', 'Longest read_post']:
        if c in df_qc.columns: df_qc[c] = df_qc[c].apply(fmt_kbp_dynamic)

    df_qc = df_qc.replace(to_replace=r'^nan[a-zA-Z\s]*$', value='', regex=True)

    def inject_hover(sample_id):
        b64_str = get_base64_image(f"03_assemblies/{sample_id}.png")
        if b64_str:
            return f'<div class="hover-container"><strong>{sample_id}</strong><img class="hover-image" src="data:image/png;base64,{b64_str}" /></div>'
        return sample_id
        
    if 'Sample' in df_qc.columns: 
        df_qc['Sample'] = df_qc['Sample'].apply(inject_hover)

    group_id = [c for c in ['Sample', 'BC'] if c in df_qc.columns]
    group_pre = [c for c in ['Median L_pre', 'Median Q_pre', 'Total reads_pre', 'Total bases_pre', 'Max Q_pre', 'Longest read_pre'] if c in df_qc.columns]
    group_post = [c for c in ['Median L_post', 'Median Q_post', 'Total reads_post', 'Total bases_post', 'Max Q_post', 'Longest read_post'] if c in df_qc.columns]
    group_ass = [c for c in ['# contigs', 'Largest contig', 'Total length', 'GC (%)', 'Depth'] if c in df_qc.columns]
    
    df_qc['|'] = ''; df_qc['| '] = ''; df_qc[' |'] = ''
    final_cols = group_id.copy()
    if group_pre: final_cols += ['|'] + group_pre
    if group_post: final_cols += ['| '] + group_post
    if group_ass: final_cols += [' |'] + group_ass
    df_qc = df_qc[final_cols]
    
    tuples = []
    for c in final_cols:
        if c in group_id: tuples.append(('Sample ID', c))
        elif c == '|': tuples.append(('sep1', ''))
        elif c in group_pre: tuples.append(('Pre-processing', c.replace('_pre', '')))
        elif c == '| ': tuples.append(('sep2', ''))
        elif c in group_post: tuples.append(('Post-processing', c.replace('_post', '')))
        elif c == ' |': tuples.append(('sep3', ''))
        elif c in group_ass: tuples.append(('Assembly', c))
    df_qc.columns = pd.MultiIndex.from_tuples(tuples)

    # ==========================================
    # TAXONOMY TABLE CONSTRUCTION
    # ==========================================
    tax_base_cols = [col for col in ['Sample', 'BC'] if not df_lista.empty and col in df_lista.columns]
    if not tax_base_cols: tax_base_cols = ['Sample']
    df_taxonomy = df_lista[tax_base_cols].copy() if not df_lista.empty and 'Sample' in df_lista.columns else pd.DataFrame(columns=tax_base_cols)
    
    if not df_tax.empty: df_taxonomy = pd.merge(df_taxonomy, df_tax, on='Sample', how='left')
    
    # Merge Kleborate (Omp mutations)
    if not df_kleb.empty and 'klebsiella_pneumo_complex__amr__Omp_mutations' in df_kleb.columns:
        df_omp = df_kleb[['Sample', 'klebsiella_pneumo_complex__amr__Omp_mutations']].copy()
        df_omp.rename(columns={'klebsiella_pneumo_complex__amr__Omp_mutations': 'Omp muts'}, inplace=True)
        # Replace the '-' that Kleborate uses so the field appears empty
        df_omp['Omp muts'] = df_omp['Omp muts'].replace('-', '')
        df_taxonomy = pd.merge(df_taxonomy, df_omp, on='Sample', how='left')

    if not df_abr.empty and 'Resistance_genes' in df_abr.columns:
        df_taxonomy = pd.merge(df_taxonomy, df_abr[['Sample', 'Resistance_genes']], on='Sample', how='left')
        df_taxonomy['Resistance_genes'] = df_taxonomy['Resistance_genes'].astype(str).str.replace(r'\s+\([\d\.]+\)', '', regex=True).replace('nan', '')

    if 'MLST' in df_taxonomy.columns:
        df_taxonomy['MLST'] = df_taxonomy['MLST'].astype(str).apply(lambda x: x[:-2] if x.endswith('.0') else x).replace('nan', '')

    for col in ['N_AMR_genes', 'AMRscore', 'VIRscore']:
        if col in df_taxonomy.columns:
            df_taxonomy[col] = df_taxonomy[col].astype(str).str.replace(',', '', regex=False).str.strip()
            df_taxonomy[col] = pd.to_numeric(df_taxonomy[col], errors='coerce').astype('Int64')
    
    df_taxonomy.drop(columns=[f'allele_{i}' for i in range(1, 8)], errors='ignore', inplace=True)
    if 'Sample' in df_taxonomy.columns: df_taxonomy['Sample'] = df_taxonomy['Sample'].apply(inject_hover)

    # Rename columns for display
    df_taxonomy = df_taxonomy.rename(columns={'Carbapenemase': 'Carba', 'ESBL': 'BLEE'}, errors='ignore')

    # Reorder columns (including Omp muts after BLEE)
    cols = list(df_taxonomy.columns)

    # Move secondary-importance columns to the end
    cols_to_end = ['Majority_genus', 'Majority_species', 'MLST_scheme', 'Possible_MLSTs', 'Possible_alleles']
    for c in cols_to_end:
        if c in cols: cols.remove(c)

    # Re-insert Omp muts after BLEE
    if 'Omp muts' in cols and 'BLEE' in cols:
        cols.remove('Omp muts')
        idx_blee = cols.index('BLEE')
        cols.insert(idx_blee + 1, 'Omp muts')
    
    # Find "Resistance_genes" to insert secondary columns after it
    if 'Resistance_genes' in cols:
        idx = cols.index('Resistance_genes')
        cols = cols[:idx+1] + [c for c in cols_to_end if c in df_taxonomy.columns] + cols[idx+1:]
    else:
        cols.extend([c for c in cols_to_end if c in df_taxonomy.columns])

    df_taxonomy = df_taxonomy[cols]

    # ==========================================
    # HTML AND CSS RENDERING
    # ==========================================
    html_template = """
    <!DOCTYPE html>
    <html lang="es">
    <head>
        <meta charset="UTF-8">
        <title>Aluminion Report</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        <link href="https://cdn.datatables.net/1.13.6/css/dataTables.bootstrap5.min.css" rel="stylesheet">
        <style>
            body {{ padding: 20px; background-color: #f8f9fa; }}
            .card {{ margin-bottom: 30px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }}
            h2 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; margin-bottom: 5px; font-weight: bold; font-size: 1.5rem; }}
            .datatable th, .datatable td {{ white-space: nowrap !important; text-align: center; vertical-align: middle; }}
            table.dataTable thead th, table.dataTable thead td, .dataframe thead th {{
                text-align: center !important; vertical-align: middle !important;
            }}
            .separator {{
                background-color: #bdc3c7 !important; min-width: 4px !important; max-width: 4px !important;
                padding: 0 !important; color: transparent !important; border: none !important;
            }}
            .hover-container {{ position: relative; display: inline-block; cursor: pointer; color: #0275d8; text-decoration: underline; }}
            .hover-image {{ 
                display: none; position: absolute; z-index: 1050; border: 3px solid #3498db; 
                border-radius: 5px; background: #fff; padding: 5px; width: 600px; 
                top: 100%; left: 50%; transform: translateX(-10%); box-shadow: 0px 10px 15px rgba(0,0,0,0.3);
            }}
            .hover-container:hover .hover-image {{ display: block; }}
            .table-responsive {{ overflow-x: auto; }}
        </style>
    </head>
    <body>
        <div class="container-fluid">
            <h1 class="mb-4" style="color: #2c3e50; font-weight: bold;">Aluminion WGS Report</h1>
            {tables_html}
        </div>
        <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.6/js/dataTables.bootstrap5.min.js"></script>
        <script>
            $(document).ready(function() {{
                $('.datatable').DataTable({{
                    "pageLength": 25, "scrollX": true,
                    "columnDefs": [{{ "orderable": false, "targets": ".separator" }}],
                    "language": {{
                        "search": "Search:", "lengthMenu": "Show _MENU_ entries per page",
                        "info": "Showing page _PAGE_ of _PAGES_",
                        "paginate": {{ "next": "Next", "previous": "Previous" }}
                    }}
                }});
            }});
        </script>
    </body>
    </html>
    """

    def df_to_html_section(df, title, item_label="samples", is_qc_table=False):
        if df.empty: return ""
        n_items = len(df)
        counter_html = f"<p style='font-weight: bold; color: #7f8c8d; margin-bottom: 15px; font-size: 1.1rem;'>N={n_items} {item_label}</p>"
        table_html = df.to_html(classes="table table-striped table-hover datatable align-middle", 
                                escape=False, index=False, justify='center', na_rep='')
        
        if is_qc_table:
            table_html = re.sub(r'(<th[^>]*)>Sample ID</th>', r'\1 style="text-align: center; font-weight: bold; color: #2c3e50; font-size: 1.1em; border-bottom: 3px solid #bdc3c7;">Sample ID</th>', table_html)
            table_html = re.sub(r'(<th[^>]*)>Pre-processing</th>', r'\1 style="text-align: center; font-weight: bold; color: #e67e22; font-size: 1.1em; border-bottom: 3px solid #e67e22;">Pre-processing</th>', table_html)
            table_html = re.sub(r'(<th[^>]*)>Post-processing</th>', r'\1 style="text-align: center; font-weight: bold; color: #27ae60; font-size: 1.1em; border-bottom: 3px solid #27ae60;">Post-processing</th>', table_html)
            table_html = re.sub(r'(<th[^>]*)>Assembly</th>', r'\1 style="text-align: center; font-weight: bold; color: #8e44ad; font-size: 1.1em; border-bottom: 3px solid #8e44ad;">Assembly</th>', table_html)
            table_html = re.sub(r'<th[^>]*>sep1</th>', '<th class="separator"></th>', table_html)
            table_html = re.sub(r'<th[^>]*>sep2</th>', '<th class="separator"></th>', table_html)
            table_html = re.sub(r'<th[^>]*>sep3</th>', '<th class="separator"></th>', table_html)  
            table_html = table_html.replace('<th></th>', '<th class="separator"></th>')
            table_html = table_html.replace('<td></td>', '<td class="separator"></td>')
        else:
            table_html = table_html.replace('<th>|</th>', '<th class="separator"></th>')
            table_html = table_html.replace('<td>|</td>', '<td class="separator"></td>')
        
        return f"""<div class="card"><div class="card-body table-responsive"><h2>{title}</h2>{counter_html}{table_html}</div></div>"""

    sections_html = ""
    sections_html += df_to_html_section(df_qc, "QC Report (Reads & Assembly)", item_label="samples", is_qc_table=True)
    sections_html += df_to_html_section(df_taxonomy, "Taxonomy & AMR Profile", item_label="samples")
    sections_html += df_to_html_section(df_copla, "Plasmids (Copla + MOBsuite)", item_label="plasmids")
    sections_html += df_to_html_section(df_integron, "Integrons (Integron Finder)", item_label="integrons")
    sections_html += df_to_html_section(df_phage, "Phages (Phastest)", item_label="phages")

    report_path = "Aluminion_Report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html_template.format(tables_html=sections_html))
    print(f"✅ HTML report generated successfully at: {report_path}")

if __name__ == "__main__":
    main()