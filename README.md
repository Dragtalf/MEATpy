# MEATpy – Microbial Ecosystem Annotation Toolkit (Python)

## About This Tool

This project, developed by Brandon Flatgard M.S., a PhD candidate under Dr. Amin Islam in the Islam Lab at Washington State University (Paul G. Allen School for Global Health, College of Veterinary Medicine), is a reimagining of the METABOLIC framework for metabolic and biogeochemical trait profiling of microbial genomes.

This project is a reimagining of the [METABOLIC](https://github.com/AnantharamanLab/METABOLIC) framework for metabolic and biogeochemical trait profiling of microbial genomes.

While inspired by METABOLIC-G and METABOLIC-C, this version was rebuilt from the ground up in Python with expanded functionality and automation, making it easier to run, integrate, and modify. Output is designed to be similar to the original METABOLIC.

---

## Key Features & Improvements

* 🔧 **Mostly automated database setup**
* 🔁 **Resumable workflows**
* ⚡ **Improved parallelization** for faster HMM searches and post-processing
* 📊 **Updated visualization support** — updated R scripts
* 🧬 **Auto-renaming of MAGs** from `nf-core/MAG`
* 🛠️ **Updated toolchains** (Diamond, CoverM, GTDB-Tk, etc.)
* 🐍 **Pure Python implementation** for maintainability and integration
* ♻️ **Smart skipping + detailed logging** for large-scale analyses

---

## Requirements

**Minimum versions tested**

```bash
python=3.13.3
HMMER=3.4
prodigal=2.6.3
sambamba=1.0.1
coverm=0.7.0
diamond=2.1.11
samtools=1.21
bowtie2=2.5.4
minimap2=2.29
gtdbtk=2.4.1
```

**Python libraries**

```bash
biopython=1.85
pandas=2.2.3
bcbio-gff=0.7.1
tqdm=4.67.1
boolean.py=5.0
openpyxl=3.1.5
requests=2.32.5
numpy=2.2.5
```

**R (for figure generation)**

```bash
r-base=4.3
r-ggplot2=3.5.2
```

---

## Install

```bash
git clone https://github.com/dragtalf/MEATpy
cd MEATpy
```

```bash
conda create -n MEATpy-env -c conda-forge \
  python=3.13.3 HMMER=3.4 prodigal=2.6.3 sambamba=1.0.1 coverm=0.7.0 \
  diamond=2.1.11 samtools=1.21 bowtie2=2.5.4 minimap2=2.29 gtdbtk=2.4.1 \
  r-base=4.3 r-ggplot2=3.5.2 \
  biopython=1.85 pandas=2.2.3 bcbio-gff=0.7.1 tqdm=4.67.1 openpyxl=3.1.5 requests=2.32.5
```

---

## KEGG Template

1. Visit: [https://www.kegg.jp/kegg-bin/get_htext?ko00002.keg](https://www.kegg.jp/kegg-bin/get_htext?ko00002.keg)
2. Click **Download json**
3. Save as `ko00002.json` inside your `templates/` directory of the MEATpy installation folder

---

## Database Setup

GTDB-Tk must be set up manually. The first run of MEATpy will configure the rest.

**Requirements**: ~140 GB for GTDB-Tk v2.4.1 reference data

Activate environment:

```bash
conda activate MEATpy-env
```

### Automatic

If you have the helper script:

```bash
download-db.sh
```

This extracts to:

```
$CONDA_PREFIX/share/gtdbtk-2.4.1/db/
```

### Manual

1. Download:

```bash
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
```

2. Extract:

```bash
tar -xvzf gtdbtk_r226_data.tar.gz -C "/path/to/target/db" --strip 1
rm gtdbtk_r226_data.tar.gz
```

3. Set env variable:

```bash
conda env config vars set GTDBTK_DATA_PATH="/path/to/target/db"
```

4. Reload env:

```bash
conda deactivate && conda activate MEATpy-env
```

---

## Command-Line Flags

| Short    | Long                    | Default        | Choices                                                         | Description                                                   |
| -------- | ----------------------- | -------------- | --------------------------------------------------------------- | ------------------------------------------------------------- |
| `-v`     | `--version`             | —              | —                                                               | Show program version                                          |
| `-t`     | `--threads`             | `20`           | —                                                               | Number of CPU threads                                         |
| `-i`     | `--input_folder`        | —              | —                                                               | Input folder containing a `fasta/` directory with genomes     |
| `-csv`   | `--input_csv`           | —              | —                                                               | CSV file with read pairs                                      |
| `-o`     | `--output_folder`       | `./MEATpy_out` | —                                                               | Output directory                                              |
| `-ko_db` | `--kofam_database_size` | `full`         | `small`, `full`                                                 | KOfam database size                                           |
| `-m`     | `--module_cutoff`       | `0.75`         | —                                                               | Cutoff value for module presence                              |
| `-pm`    | `--prodigal_method`     | `meta`         | `meta`, `single`                                                | Prodigal mode for ORF calling                                 |
| `-rt`    | `--reads_type`          | `metaG`        | `metaG`, `metaT`                                                | Type of omic reads                                            |
| `-st`    | `--sequencing_type`     | `illumina`     | `illumina`, `pacbio`, `pacbio_hifi`, `pacbio_asm20`, `nanopore` | Sequencing type                                               |
| `-tax`   | `--taxonomy`            | `phylum`       | `phylum`, `class`, `order`, `family`, `genus`, `species`, `bin` | Taxonomic level for MW-score                                  |
| `-c`     | `--correction`          | `GTDB`         | `GTDB`, `NCBI`                                                  | Taxonomy correction type                                      |
| —        | `--skip_to_diagrams`    | —              | —                                                               | Start directly at diagram step                                |
| —        | `--rename`              | —              | —                                                               | Rename MAGs from `nf-core/mag` using GTDB/NCBI classification |
| —        | `--revert_names`        | —              | —                                                               | Revert renamed MAGs to originals                              |
| —        | `--oops`                | —              | —                                                               | Fix if you forgot `--rename`                                  |

---

## Input
1. MAGs or samples in .fasta format in a "fasta" folder. The program will create the other relavent files.
    MAGs prepared using nf-core/MAG can be copied to the fasta folder and can be optionally renamed using the program commands.
    example: ~/INPUT_FOLDER/fasta/ALL_YOUR_FASTA_FILES
2. A .csv containing the path to your reads
    example: ~/reads/Sample_R1.fastq.gz, ~/reads/Sample_R2.fastq.gz

## Running

Activate environment:

```bash
conda activate MEATpy-env
```

**Example:**

```bash
python start.py -i ~/MEATpy_input/ -csv ~/sample_reads.csv -o ~/MEATpy_Out -tax genus -c NCBI -t 20
```

---

## Recommended Citation

If you use this software, please cite **both** this tool and the original METABOLIC paper.

**MEATpy – Microbial Ecosystem Annotation Toolkit (Python)**
## Citation

> Flatgard, B. M. & Islam, M. A. (2025). *MEATpy – Microbial Ecosystem Annotation Toolkit (Python)*. [doi.org/10.5281/zenodo.17317008](https://doi.org/10.5281/zenodo.17317008).  
> GitHub: [https://github.com/Dragtalf/MEATpy](https://github.com/Dragtalf/MEATpy)

**METABOLIC**
> Zhou, Z., Tran, P.Q., Breister, A.M., *et al.* (2022).
> METABOLIC: high-throughput profiling of microbial genomes for functional traits, metabolism, biogeochemistry, and community-scale functional networks.
> *Microbiome* **10**, 33. [doi.org/10.1186/s40168-021-01213-8](https://doi.org/10.1186/s40168-021-01213-8)

---

## License & Credits

This project is released under the [GPL-3.0 License](LICENSE).
Portions of data formats and logic are adapted from METABOLIC (Zhou *et al.*, 2022).
