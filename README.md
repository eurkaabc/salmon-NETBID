# Salmon-NETBID2 Workflow

## Introduction

**Salmon** is a fast and efficient tool for transcript quantification from RNA-seq data.  
It requires a set of target transcripts (either from a reference genome annotation or from a de-novo assembly) to quantify.  

To run **Salmon**, you will need:  
- A FASTA file containing the reference transcripts  
- One or more FASTA/FASTQ files containing the RNA-seq reads  

Optionally, Salmon can also make use of pre-computed alignments (SAM/BAM files) to the transcripts instead of raw reads.

In this workflow, we assume that the **raw FASTQ files** are stored in:  
```bash
/mnt/sda/Public/Project/collabration/AoLab/20250821/rawdata
```

---

<img width="658" height="428" alt="image" src="https://github.com/user-attachments/assets/9e9eaea0-eac8-4cab-96e9-2a3e18d32999" />

## Installation
## Step 0. 整理FASTQ 
```bash
#
mkdir -p "/mnt/sda/Public/Project/collabration/AoLab/20250821/1.data"
find "/mnt/sda/Public/Project/collabration/AoLab/20250821/rawdata" -type f -name "*.fastq.gz" -exec cp -n {} "/mnt/sda/Public/Project/collabration/AoLab/20250821/1.data/" \;
```

## 🧪 Step 1. Build Salmon Index (Mouse GRCm38, Gencode vM23)
```bash
cd /mnt/sda/Public/Database/salmon_usage

# 1. Download transcriptome and genome files
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz

# 2. Install Salmon (make sure you are inside a conda environment)
conda install --channel bioconda salmon

# 3. Generate the decoy list (chromosome IDs from the genome fasta)
grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

# 4. Build the gentrome file (transcriptome + genome)
cat gencode.vM23.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome.fa.gz

# 5. Build the Salmon index
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_M38_index --gencode
```
## Step 1 (Alternative). Build Salmon Index (Mouse GRCm39, Ensembl Release 115)
```bash
cd /mnt/sda/Public/Database/salmon_usage

# 1. Download transcriptome (cDNA) and genome
wget ftp://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# 2. Define file paths and output directory
CDNA="/mnt/sda/Public/Database/salmon_usage/Mus_musculus.GRCm39.cdna.all.fa.gz"
GENOME="/mnt/sda/Public/Database/salmon_usage/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
INDEX_DIR="/home/Weixin.Zhang/rnaseq_index/salmon_M39_index"
THREADS=16

mkdir -p "$(dirname "$INDEX_DIR")"
cd "$(dirname "$INDEX_DIR")"

# 3. Generate decoy list
zgrep -h '^>' "$GENOME" | cut -d ' ' -f1 | sed 's/^>//' > decoys.txt

# 4. Build gentrome (cDNA + genome)
cat <(zcat "$CDNA") <(zcat "$GENOME") | gzip -c > gentrome.fa.gz

# 5. Build the Salmon index
salmon index -t gentrome.fa.gz -d decoys.txt -i "$INDEX_DIR" -k 31 -p "$THREADS" > salmon_index.log 2>&1
```


## 🧪 Step 2. Salmon Quantification (Batch, Paired-End)

Quantify all paired-end FASTQ samples under your project.  
This script auto-detects library type (`-l A`), enables recommended bias corrections, and writes a per-sample summary.

> **Prerequisites**
> - Salmon installed (see Step 1)
> - `jq` for parsing `meta_info.json`  
>   - Install via conda: `conda install -c conda-forge jq`
> - `column` (usually in `bsdmainutils`/`util-linux`, available by default on most Linux distros)

```bash
# --- Paths & Settings ---
PROJECT_ROOT="/mnt/sda/Public/Project/collabration/AoLab/20250821"
RAW="$PROJECT_ROOT/1.data"
INDEX="/mnt/sda/Public/Database/salmon_usage/salmon_M38_index"        # e.g., salmon_M38_index or salmon_M39_index
OUTROOT="$PROJECT_ROOT/2.salmon"

THREADS=12

mkdir -p "$OUTROOT"

# --- Quantify all paired-end samples ---
for R1 in "$RAW"/*_R1.fastq.gz; do
  S=$(basename "$R1" _R1.fastq.gz)
  R2="$RAW/${S}_R2.fastq.gz"

  if [ ! -f "$R2" ]; then
    echo "[WARN] Missing mate for sample '$S': $R2 (skipping)"
    continue
  fi

  OUT="$OUTROOT/$S"
  mkdir -p "$OUT"

  echo "[INFO] Quantifying sample: $S"
  salmon quant \
    -i "$INDEX" \
    -l A \
    -1 "$R1" \
    -2 "$R2" \
    -p "$THREADS" \
    --validateMappings \
    --gcBias \
    --seqBias \
    -o "$OUT"
done
```

**Parameters Explained**

| Parameter            | Description                                                                                        |
| -------------------- | -------------------------------------------------------------------------------------------------- |
| `-i "$INDEX"`        | Path to Salmon index (from Step 1).                                                                |
| `-l A`               | Automatically detect library type (works for most datasets).                                       |
| `-1 "$R1" -2 "$R2"`  | Paired-end FASTQ files (read 1 and read 2).                                                        |
| `-p "$THREADS"`      | Number of CPU threads to use (adjust to your system).                                              |
| `--validateMappings` | More accurate mapping algorithm (quasi-mapping + validation). Recommended for decoy-aware indexes. |
| `--gcBias`           | Corrects for GC-content bias in RNA-seq.                                                           |
| `--seqBias`          | Corrects for random hexamer priming bias in RNA-seq.                                               |
| `-o "$OUT"`          | Output directory for the sample. Contains `quant.sf` (expression estimates).                       |


```bash
# --- Build a small summary table (one row per sample) ---
SUMMARY_TSV="$OUTROOT/percent_mapped_summary.tsv"
TMP_SUMMARY="$OUTROOT/.percent_mapped_summary.tmp"

# header
echo -e "sample\tlibrary_type\tpercent_mapped\tnum_mapped\tnum_decoy_fragments" > "$TMP_SUMMARY"

# rows
for m in "$OUTROOT"/*/aux_info/meta_info.json; do
  S=$(basename "$(dirname "$(dirname "$m")")")
  jq -r "[\"${S}\", .library_type, .percent_mapped, .num_mapped, .num_decoy_fragments] | @tsv" "$m" \
    >> "$TMP_SUMMARY"
done

mv "$TMP_SUMMARY" "$SUMMARY_TSV"

# pretty-print to terminal (optional)
echo
echo "[INFO] Summary written to: $SUMMARY_TSV"
echo "[INFO] Preview:"
column -t -s $'\t' "$SUMMARY_TSV"
```

<img width="6256" height="4167" alt="image" src="https://github.com/user-attachments/assets/fea832e0-f2b2-42e2-a966-9d0be8072c86" />
