# nanoesst

A bioinformatics tool for Nanopore sequencing processing, mapping 3rd-generation reads using `fastplong`, determining species composition with `sylph`, mapping to ESKAPEE reference genomes via `minimap2`, and identifying ST types using `pymlst`.

## Installation & Environment

Create a conda environment with all dependencies, including `kraken2` and `bracken` for the new profiling features:

```bash
conda create -n nanoesst fastplong sylph minimap2 pymlst samtools pigz kraken2 bracken
conda activate nanoesst
```

Then, install this package directly from the source code:

```bash
cd nanoesst-main
pip install -e .
```

## Usage

The tool requires you to choose an identification algorithm using `-a sylph` or `-a kraken`. Based on your choice, you must provide the corresponding database path (`-syldb` or `-krakendb`). 

### 1. Process Mode (Single Sample)

**Using Sylph:**
```bash
nanoesst process -i barcode01.fastq.gz -n SK-1 -a sylph -syldb path/to/database.syldb -t 16
```

**Using Kraken2:**
```bash
nanoesst process -i barcode01.fastq.gz -n SK-1 -a kraken -krakendb path/to/kraken_db -t 16
```

### 2. Batch Mode (Multiple Samples)

**Using Kraken2:**
```bash
nanoesst batch -i ./raw_fastq_dir/ -n mapping.txt -a kraken -krakendb path/to/kraken_db -t 16
```

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
