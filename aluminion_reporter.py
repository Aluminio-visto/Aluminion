#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import base64
import os
import sys

def get_base64_image(image_path):
    """Lee una imagen PNG y la convierte a una cadena Base64."""
    if os.path.exists(image_path):
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    return None

def load_and_standardize(file_path, sep=',', key_col=None):
    """Carga un archivo (CSV, TSV o Excel), estandariza la columna llave a 'Sample' y lo devuelve."""
    if not os.path.exists(file_path):
        print(f"Aviso: No se encontró el archivo {file_path}")
        return pd.DataFrame()
    
    try:
        if file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path)
        else:
            df = pd.read_csv(file_path, sep=sep)
        
        # Renombrar la columna llave si se especifica
        if key_col and key_col in df.columns:
            df.rename(columns={key_col: 'Sample'}, inplace=True)
            
        # Asegurar que 'Sample' sea texto puro para evitar fallos al cruzar datos
        if 'Sample' in df.columns:
            df['Sample'] = df['Sample'].astype(str).str.strip()
            
        return df
    except Exception as e:
        print(f"Error al leer {file_path}: {e}")
        return pd.DataFrame()

def main():
    # 1. Definir directorio de trabajo (por defecto el actual)
    work_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    os.chdir(work_dir)
    print(f"Generando informe en: {os.getcwd()}")

    # 2. Cargar y estandarizar todas las tablas
    df_lista    = load_and_standardize("lista_seq.tsv", sep='\t', key_col='ID')
    df_qcreads  = load_and_standardize("QC_reads.csv", sep='\t', key_col='Sample')
    df_qcass    = load_and_standardize("QC_assembly.csv", sep='\t', key_col='Samples')
    df_tax      = load_and_standardize("taxonomy.xlsx", key_col='Sample')
    
    df_copla    = load_and_standardize("copla_modif.csv", sep=',', key_col='Sample')
    df_integron = load_and_standardize("integron_summary.csv", sep=',', key_col='Sample')
    df_phage    = load_and_standardize("phage_summary.csv", sep=',', key_col='sample')

    # 3. Construir la Tabla Maestra (Lista + QC + Taxonomía)
    master_df = df_lista.copy()
    if not df_qcreads.empty:
        master_df = pd.merge(master_df, df_qcreads, on='Sample', how='left')
    if not df_qcass.empty:
        master_df = pd.merge(master_df, df_qcass, on='Sample', how='left')
    if not df_tax.empty:
        master_df = pd.merge(master_df, df_tax, on='Sample', how='left')

    # 4. Inyectar el efecto Hover con las imágenes de Bandage en la columna Sample
    if 'Sample' in master_df.columns:
        def inject_hover(sample_id):
            img_path = f"03_assemblies/{sample_id}.png"
            b64_str = get_base64_image(img_path)
            if b64_str:
                # El CSS se encarga de mostrar .hover-image solo al pasar el ratón
                return f'<div class="hover-container"><strong>{sample_id}</strong><img class="hover-image" src="data:image/png;base64,{b64_str}" /></div>'
            return sample_id
        
        master_df['Sample'] = master_df['Sample'].apply(inject_hover)

    # 5. Configurar el renderizado HTML
    # Plantilla base con Bootstrap 5 y DataTables
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
            h2 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; margin-bottom: 20px; }}
            
            /* CSS para el hover de las imágenes de Bandage */
            .hover-container {{ position: relative; display: inline-block; cursor: pointer; color: #0275d8; text-decoration: underline; }}
            .hover-image {{ 
                display: none; 
                position: absolute; 
                z-index: 1050; 
                border: 3px solid #3498db; 
                border-radius: 5px;
                background: #fff; 
                padding: 5px; 
                width: 600px; /* Tamaño de la imagen flotante */
                top: 100%; 
                left: 50%; 
                transform: translateX(-10%);
                box-shadow: 0px 10px 15px rgba(0,0,0,0.3);
            }}
            .hover-container:hover .hover-image {{ display: block; }}
            
            /* Ajustes para tablas largas */
            .table-responsive {{ overflow-x: auto; }}
            th {{ white-space: nowrap; }}
        </style>
    </head>
    <body>
        <div class="container-fluid">
            <h1 class="mb-4">🧬 Aluminion WGS Report</h1>
            
            {tables_html}
            
        </div>
        
        <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.6/js/dataTables.bootstrap5.min.js"></script>
        <script>
            $(document).ready(function() {{
                $('.datatable').DataTable({{
                    "pageLength": 25,
                    "scrollX": true,
                    "language": {{
                        "search": "Buscar:",
                        "lengthMenu": "Mostrar _MENU_ registros por página",
                        "info": "Mostrando página _PAGE_ de _PAGES_",
                        "paginate": {{ "next": "Siguiente", "previous": "Anterior" }}
                    }}
                }});
            }});
        </script>
    </body>
    </html>
    """

    # Función para convertir dataframe a HTML de Bootstrap
    def df_to_html_section(df, title):
        if df.empty:
            return ""
        table_html = df.to_html(classes="table table-striped table-hover datatable", escape=False, index=False, justify='center')
        return f"""
        <div class="card">
            <div class="card-body table-responsive">
                <h2>{title}</h2>
                {table_html}
            </div>
        </div>
        """

    # 6. Ensamblar todas las secciones
    sections_html = ""
    sections_html += df_to_html_section(master_df, "📊 Master Table (Taxonomy & QC)")
    sections_html += df_to_html_section(df_copla, "🧬 Plasmids (Copla)")
    sections_html += df_to_html_section(df_integron, "🧩 Integrons (Integron Finder)")
    sections_html += df_to_html_section(df_phage, "🦠 Phages (Phastest)")

    # 7. Guardar el archivo final
    report_path = "Aluminion_Report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html_template.format(tables_html=sections_html))
    
    print(f"✅ Informe HTML generado con éxito en: {report_path}")

if __name__ == "__main__":
    main()
