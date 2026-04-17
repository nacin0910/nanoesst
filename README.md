# nanoesst

A bioinformatics tool for Nanopore sequencing processing, mapping 3rd-generation reads using `fastplong`, determining species composition with `sylph`, mapping to ESKAPEE reference genomes via `minimap2`, and identifying ST types using `pymlst`.

## Installation & Environment

First, set up the conda environment with the required dependencies (note: `samtools` is also required for BAM manipulation):

```bash
conda create -n nanoesst fastplong sylph minimap2 pymlst samtools pigz
conda activate nanoesst
```

Then, install this package directly from the source code:

```bash
cd nanoesst-main
pip install -e .
```

> **Important**: Run the tool from the `nanoesst-main` directory, as it expects the reference genomes and MLST databases to be located in the local `db/` folder.

## Usage

The tool supports two modes: `process` (for a single file) and `batch` (for a folder of files).

### 1. Process Mode (Single Sample)
Use this mode to process a single fastq.gz file.

```bash
nanoesst process -i barcode01.fastq.gz -n SK-1 -db path/to/database.syldb -t 16
```
* `-i`: Path to the input fastq.gz file.
* `-n`: Sample name (will be used as the prefix for all output files).
* `-db`: Path to the sylph database.

### 2. Batch Mode (Multiple Samples)
Use this mode to process a directory containing multiple fastq.gz files.

```bash
nanoesst batch -i ./raw_fastq_dir/ -n mapping.txt -db path/to/database.syldb -t 16
```
* `-i`: Directory containing all `.fastq.gz` files.
* `-n`: A two-column mapping file (tab or space separated). The first column is the barcode (e.g., `barcode01`), and the second column is the sample name (e.g., `SK-1`).

**Example `mapping.txt`:**
```text
barcode01   SK-1
barcode02   SK-2
barcode03   SK-3
```

## Pipeline Workflow

1. **Naming**: Assigns the user-provided sample name.
2. **QC**: Runs `fastplong` to output `1.fastplong/$sample_clean.fastq.gz`.
3. **Profiling**: Runs `sylph profile` to output `2.sylph_out/$sample_profile.tsv`. It then parses the `Contig_name` and `Taxonomic_abundance` columns to check if any ESKAPEE pathogens are present (>1% abundance).
4. **Mapping**: If present, maps clean reads to the pathogen reference genome using `minimap2` and extracts them using `samtools` to `3.minimap2/$sample_pathogen.fastq.gz`.
5. **Typing**: Runs `claMLST` on the mapped reads and writes ST results to `4.pymlst/$sample_pathogen_result.txt`.
