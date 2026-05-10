# Aluminion 🧬

**Aluminion** is an automated, highly-scalable, and modular Bioinformatics pipeline for the assembly, annotation, and comprehensive characterization of bacterial whole-genome sequencing (WGS) data from Oxford Nanopore Technologies (MinION).

## 🚀 Key Features

* **End-to-End Processing**: From raw Nanopore `fastq_pass` reads to polished assemblies.
* **Modular Architecture**: Run the full suite or skip specific heavy modules (`--skip-phages`, `--skip-plasmids`, etc.) to save time.
* **Rich Annotation Profile**: Integrated detection of AMR genes, Plasmids, Phages, Integrons, and MLST schemes.
* **Interactive HTML Reports**: Automatic generation of interactive DataTables with hoverable assembly graphs (Bandage).

## 🛠️ Tools Included

Aluminion chains together industry-standard bioinformatics software:
* **QC & Filtering**: NanoPlot, Chopper
* >> Assembly: Flye, Circlator, Bandage, Quast
* **Taxonomy & Typing**: Kraken2, GAMBIT, Kleborate, ECTyper, MLST
* **Annotation & Genomics**: Bakta, Abricate, MOB-suite, Copla, Integron_Finder, Phastest (Docker), ISfinder

## 📥 Installation

Aluminion relies heavily on Conda environments to avoid dependency conflicts. 

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/YourOrganization/aluminion.git](https://github.com/YourOrganization/aluminion.git)
   cd aluminion
   chmod +x aluminion.sh

```

2. **Set up Conda Environments:**
Use the provided `environment.yml` to set up the base environment.
```bash
conda env create -f envs/environment.yml

```


*(Note: Sub-environments for Copla, Integron_Finder, and Kleborate are handled by the pipeline).*
3. **Database setup:**
You will need to download and configure the required databases for Kraken2, GAMBIT, Bakta, Copla, and ISfinder in a dedicated directory.

## ⚙️ Usage

Run the pipeline using the orchestrator bash script. A typical execution looks like this:

```bash
./aluminion.sh \\
  -r RUN_FOLDER_NAME \\
  -b /path/to/Databases \\
  -t 30 \\
  -l /path/to/lista_seq.tsv

```

### 🚦 Skipping Modules (Time-saving Flags)

You can skip computationally heavy or unneeded modules by passing flags:

* `--skip-phages`: Skips Phastest execution and parsing.
* `--skip-integrons`: Skips Integron_Finder execution and parsing.
* `--skip-plasmids`: Skips Copla plasmid clustering.
* `--skip-typing`: Skips MLST, Kleborate, and ECTyper consolidation.
* `--skip-kraken`: Skips Kraken2 taxonomy parsing.
* `--skip-abr`: Skips Abricate AMR parsing.

Example:

```bash
./aluminion.sh -r BAC_2025 -b /mnt/db -t 30 --skip-phages --skip-integrons

```

## 📊 Outputs

Upon completion, Aluminion generates a structured working directory:

* `01_reads/` & `02_filter/`: Cleaned and filtered reads + NanoPlot QC.
* `03_assemblies/`: Flye assemblies, Bandage plots, and Quast reports.
* `04_taxonomies/`: Kraken2, GAMBIT, and Kleborate logs.
* `08_Anotacion/`: Bakta annotations, MOB-suite extracts, and Abricate tables.
* `09_phages/` `05_plasmids/` & `11_integrons/`: Specific mobile element outputs.
* **`Aluminion_Report.html`**: The final interactive report summarizing QC, Taxonomy, AMR, Plasmids, Phages, and Integrons.

## 🤝 Contributing

Contributions, issues, and feature requests are welcome! We plan to expand the pipeline to include Virulence Factors and Anti-Phage defense systems in the near future.
"""



```