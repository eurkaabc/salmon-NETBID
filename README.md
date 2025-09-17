# Salmon-NETBID2 Workflow

## 📖 Introduction

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

## ⚙️ Installation

Make sure you have conda installed, then install Salmon via Bioconda:
```bash
conda install --channel bioconda salmon

```
## Table of Contents
- [🧬 Step 1. Build Salmon Index](#-step-1-build-salmon-index)
- [🔬 Step 2. Salmon Quantification (Batch, Paired-End)](#-step-2-salmon-quantification-batch-paired-end)
- [🗂️ Step 3. Prepare Files for R Analysis](#-step-3-prepare-files-for-r-analysis)
  - [🧩 Build tx2gene.csv mapping file](#-build-tx2genecsv-mapping-file)
- [🧠 Step 4. Load in R & Gene ID Mapping](#-step-4-load-in-r--gene-id-mapping)
- [🕸️ Step 5. Run SJARACNe for Network Inference (Server)](#-step-5-run-sjaracne-for-network-inference-server)
- [🧭 Step 6. NetBID2 Hidden Driver Estimation](#-step-6-netbid2-hidden-driver-estimation)
- [🧮 Step 7. Differential Expression/Activity (KO vs WT) & Master Table](#-step-7-differential-expressionactivity-ko-vs-wt--master-table)
- [🚀 Step 8. Advanced Analysis (Volcano, GSEA, Enrichment, Heatmaps & Network)](#-step-8-advanced-analysis-volcano-gsea-enrichment-heatmaps--network)

## 📂 Step 0. Organize FASTQ Files
```bash
mkdir -p "/mnt/sda/Public/Project/collabration/AoLab/20250821/1.data"
find "/mnt/sda/Public/Project/collabration/AoLab/20250821/rawdata" -type f -name "*.fastq.gz" -exec cp -n {} "/mnt/sda/Public/Project/collabration/AoLab/20250821/1.data/" \;
```

## 🧬 Step 1. Build Salmon Index
**Choose the index according to your species and annotation source:**
<details> <summary><strong>Option A: Mouse GRCm38 (GENCODE vM23)</strong></summary>

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
</details>




<details> <summary><strong>Option B: Mouse GRCm39 (Ensembl Release 115)</strong></summary>
  
```bash
cd /mnt/sda/Public/Database/salmon_usage

# 1. Download transcriptome (cDNA) and genome
wget ftp://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# 2. Define paths
CDNA="/mnt/sda/Public/Database/salmon_usage/Mus_musculus.GRCm39.cdna.all.fa.gz"
GENOME="/mnt/sda/Public/Database/salmon_usage/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
INDEX_DIR="/home/Weixin.Zhang/rnaseq_index/salmon_M39_index"
THREADS=16

mkdir -p "$(dirname "$INDEX_DIR")"
cd "$(dirname "$INDEX_DIR")"

# 3. Generate decoy list
zgrep -h '^>' "$GENOME" | cut -d ' ' -f1 | sed 's/^>//' > decoys.txt

# 4. Build gentrome 
cat <(zcat "$CDNA") <(zcat "$GENOME") | gzip -c > gentrome.fa.gz

# 5. Build the Salmon index
salmon index -t gentrome.fa.gz -d decoys.txt -i "$INDEX_DIR" -k 31 -p "$THREADS" > salmon_index.log 2>&1
```
</details>


<details> <summary><strong>Option C: Human GENCODE v45</strong></summary>

Note: Build a Salmon index using GENCODE v45 human transcriptome (no decoy-aware mode in this example).
```bash
cd /mnt/sda/Public/Database/salmon_usage
mkdir human_gencode_v45
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.transcripts.fa.gz

salmon index \
  -i human_gencode_v45 \
  -t gencode.v45.transcripts.fa.gz \
  --gencode \
  -k 31 \
  -p 8
```
</details>





## 🔬 Step 2. Salmon Quantification (Batch, Paired-End)

This step runs Salmon quantification across all paired-end FASTQ samples, producing per-sample `quant.sf` expression files and a summary table of mapping statistics.

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
**✅ Expected Outputs**

After running this step, you should have:

| File / Folder                              | Description                                         |
| ------------------------------------------ | --------------------------------------------------- |
| `2.salmon/<sample>/quant.sf`               | Transcript-level quantification results             |
| `2.salmon/<sample>/aux_info/meta_info.json`| Mapping statistics (used for summary)               |
| `2.salmon/percent_mapped_summary.tsv`      | Summary table (sample name, mapping %, counts, etc.)|



## 🗂️ Step 3. Prepare Files for R Analysis

Once all samples have been quantified with Salmon, we need to organize the results into a format that can be easily imported into **R** (e.g., with `tximport` for NETBID2).

## 🔧 Prepare quantification outputs
```bash
# Define project paths
RDIR="$PROJECT_ROOT/3.R_analysis"
mkdir -p "$RDIR"

QUANTS_DIR="${PROJECT_ROOT}/2.salmon"

# Move into the Salmon quantification directory
cd "$QUANTS_DIR"

# 1. Compress all quant.sf files into one archive (for backup or sharing)
zip "$RDIR/quant.sf.zip" $(find . -name quant.sf)

# 2. List all sample folder names that contain quant.sf
for q in */quant.sf; do
  basename "$(dirname "$q")"
done | sort > "$RDIR/sampleFile"

# 3. Generate salmon.output (mapping sample → quant.sf path)
awk -v QUANTS_DIR="$QUANTS_DIR" '{printf "%s\t%s/%s/quant.sf\n",$1,QUANTS_DIR,$1}' "$RDIR/sampleFile" > "$RDIR/salmon.output"

echo "[INFO] R analysis files written to: $RDIR"
ls -lh "$RDIR"

# 4. (Optional) unzip the archive inside the R analysis directory
cd "$RDIR"
unzip -o quant.sf.zip
```



**🧩Build tx2gene.csv mapping file**

Depending on whether your quant.sf transcript IDs retain version numbers (e.g., ENSMUST00000193812.1) or not, generate a matching tx2gene.csv.


## A) Keep version (if `quant.sf` IDs like `ENSMUST00000193812.1`):
```bash
CDNA="/mnt/sda/Public/Database/salmon_usage/gencode.vM23.transcripts.fa.gz"
RDIR="$PROJECT_ROOT/3.R_analysis"

zgrep '^>' "$CDNA" \
| awk -F'|' 'BEGIN{OFS=","; print "transcript","gene"}{
  tx=$1; g=$2; sub(/^>/,"",tx);   # keep version
  sub(/^>/,"",g);                 # keep gene version
  print tx,g
}' > "$RDIR/tx2gene.csv"

head "$RDIR/tx2gene.csv"

```


## B) Remove version (if `quant.sf` lacks them):
```bash

CDNA="/mnt/sda/Public/Database/salmon_usage/gencode.vM23.transcripts.fa.gz"
RDIR="$PROJECT_ROOT/3.R_analysis"

zgrep '^>' "$CDNA" \
| awk -F'|' 'BEGIN{OFS=","; print "transcript","gene"}{
  tx=$1; g=$2; sub(/^>/,"",tx);         # ENSMUST... .xx
  sub(/\..*$/,"",tx); sub(/\..*$/,"",g) # strip version
  print tx,g
}' > "$RDIR/tx2gene.csv"

head "$RDIR/tx2gene.csv"

```


**✅ Expected Outputs**

After this step, you should have:
| File / Folder                | Description                                                                 |
| ---------------------------- | --------------------------------------------------------------------------- |
| `3.R_analysis/quant.sf.zip`  | Archive of all `quant.sf` files (backup/sharing).                           |
| `3.R_analysis/sampleFile`    | List of sample names (one per line).                                        |
| `3.R_analysis/salmon.output` | Mapping of `sample → quant.sf` file path (used by `tximport` / NetBID2).    |
| `3.R_analysis/tx2gene.csv`   | Transcript-to-gene mapping (with or without version numbers, as generated). |









## 🧠 Step 4. Load in R & Gene ID Mapping

### Stage I: Load & Mapping
```bash
# ---- Load necessary libraries ----
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(NetBID2); library(Biobase); library(edgeR)
})

# ---- Define Paths ----
project_dir <- "/mnt/sda/Public/Project/collabration/AoLab/20250821"  # **Modify this path if needed**
salmon_dir  <- file.path(project_dir, "2.salmon")
rdir        <- file.path(project_dir, "3.R_analysis")

# ---- Set working directory ----
setwd(rdir)  # **Set to R analysis directory**

# ---- Create directories if they do not exist ----
dir.create(rdir, showWarnings = FALSE, recursive = TRUE)

# ---- Sample table (already prepared) ----
sampleFile <- "sampleFile"
samples <- readLines(sampleFile)

# ---- Template metadata (fill WT/KO later) ----
coldata <- data.frame(
  sample = samples,
  group  = NA_character_,   # **Fill with WT/KO (or your groups)**
  stringsAsFactors = FALSE
)
write.csv(coldata, "sample_metadata.template.csv", row.names = FALSE)

# ---- Load completed metadata ----
meta <- read.csv(file.path(rdir, "sample_metadata.csv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$sample

# ---- tx2gene (match your quant.sf: keep version if quant.sf has version) ----
tx2 <- read.csv(file.path(rdir, "tx2gene.csv"),
                header = FALSE, stringsAsFactors = FALSE)
names(tx2) <- c("TXNAME","GENEID")


# ---- 1) Load gene-level ExpressionSet from Salmon ----
eSet_expGene <- load.exp.RNASeq.demoSalmon(
  salmon_dir         = salmon_dir,
  tx2gene            = tx2,           # or a file path string
  use_phenotype_info = meta,
  use_sample_col     = "sample",
  use_design_col     = "group",
  merge_level        = "gene",
  return_type        = "eset"         # **case-sensitive**
)
saveRDS(eSet_expGene, file.path(rdir, "Gene.eset.rawcounts.rds"))
cat("[OK] eSet维度: ", paste(dim(eSet_expGene), collapse=" x "), "\n", sep="")

# ---- 2) Map Ensembl Gene IDs → Gene Symbols ----
# Choose species dataset: set to mouse or human explicitly
species_dataset <- "mmusculus_gene_ensembl"     # **use "hsapiens_gene_ensembl" for human**

transfer_tab <- get_IDtransfer2symbol2type(
  from_type      = "ensembl_gene_id",
  dataset        = species_dataset,
  use_level      = "gene",
  ignore_version = TRUE
)

# Prefer species-specific symbol column; fallback to external_gene_name
to_col <- if ("mgi_symbol" %in% colnames(transfer_tab)) {
  "mgi_symbol"
} else if ("hgnc_symbol" %in% colnames(transfer_tab)) {
  "hgnc_symbol"
} else {
  "external_gene_name"
}

## A) strip version from featureNames and ensure uniqueness
old_fn   <- featureNames(eSet_expGene)
fn_nov   <- sub("\\.\\d+$","", old_fn)
fn_nov_u <- make.unique(fn_nov, sep = ".dup")
featureNames(eSet_expGene) <- fn_nov_u

## B) sync featureData (store Ensembl IDs without version)
fd <- data.frame(ensembl_gene_id = fn_nov,
                 row.names = fn_nov_u,
                 stringsAsFactors = FALSE)
featureData(eSet_expGene) <- new("AnnotatedDataFrame", data = fd)

## C) apply mapping & merge duplicates by median
eSet_expGene2 <- update_eset.feature(
  use_eset         = eSet_expGene,
  use_feature_info = transfer_tab,
  from_feature     = "ensembl_gene_id",
  to_feature       = to_col,
  merge_method     = "median"
)

saveRDS(eSet_expGene2, file.path(rdir, "Gene.eset.mapped.rds"))
cat("[OK] After mapping dims: ", paste(dim(eSet_expGene2), collapse=" x "), "\n", sep="")
```


### Stage II: QC + Save
```bash
# ---- Init parameter container ----
network.par <- list()

# ---- QC output directory ----
network.par$out.dir.QC <- file.path(rdir, "QC")
dir.create(network.par$out.dir.QC, showWarnings = FALSE, recursive = TRUE)

# ---- Put mapped eSet into network.par ----
network.par$net.eset <- eSet_expGene2

# ---- Strictly align metadata order to expression columns ----
stopifnot(all(colnames(exprs(network.par$net.eset)) %in% meta$sample))
meta <- meta[colnames(exprs(network.par$net.eset)), , drop = FALSE]
stopifnot(all(colnames(exprs(network.par$net.eset)) == rownames(meta)))

# ---- QC (intgroup must exist in metadata; here we use "group") ----
draw.eset.QC(
  network.par$net.eset,
  outdir          = network.par$out.dir.QC,
  intgroup        = "group",
  do.logtransform = FALSE,
  prefix          = "beforeQC_",
  generate_html   = FALSE
)

# ---- Define data save directory (use rdir) ----
network.par$out.dir.DATA <- rdir
dir.create(network.par$out.dir.DATA, showWarnings = FALSE, recursive = TRUE)

# ---- Save current step ----
NetBID.saveRData(network.par = network.par, step = "exp-load")

```

### Stage III: Normalization

```bash
# ---- Extract expression matrix ----
mat <- exprs(network.par$net.eset)

# ---- Low-expression filter (≤ 5th percentile in ≥90% samples) ----
choose1 <- apply(mat <= quantile(mat, probs = 0.05), 1, sum) <= ncol(mat) * 0.90
cat("[FILTER] Low-expression kept/dropped:\n"); print(table(choose1))
mat <- mat[choose1, , drop = FALSE]

# ---- Update ExpressionSet ----
net_eset <- generate.eset(
  exp_mat         = mat,
  phenotype_info  = pData(network.par$net.eset)[colnames(mat), , drop = FALSE],
  feature_info    = fData(network.par$net.eset)[rownames(mat), , drop = FALSE],
  annotation_info = annotation(network.par$net.eset)
)

# ---- Update network.par ----
network.par$net.eset <- net_eset

# ---- QC after normalization ----
stopifnot("group" %in% colnames(pData(network.par$net.eset)))  # ensure 'group' exists in metadata
draw.eset.QC(
  network.par$net.eset,
  outdir          = network.par$out.dir.QC,
  intgroup        = "group",    # ← must match your metadata column
  do.logtransform = FALSE,
  prefix          = "afterQC_",
  generate_html   = FALSE
)

# ---- Save current step ----
NetBID.saveRData(network.par = network.par, step = "exp-QC")

# ---------------- Optional: Sample clustering check ----------------

intgroup <- "group"  # must match a column in your metadata
mat <- exprs(network.par$net.eset)

# ---- High-variance gene filter via IQR ----
choose1 <- IQR.filter(exp_mat = mat, use_genes = rownames(mat), thre = 0.8)
cat("[FILTER] High-variance kept/dropped:\n"); print(table(choose1))
mat <- mat[choose1, , drop = FALSE]

# ---- K-means clustering vs observed labels ----
pred_label <- draw.emb.kmeans(
  mat        = mat,
  all_k      = NULL,
  obs_label  = get_obs_label(pData(network.par$net.eset), intgroup),
  pre_define = c('WNT'='blue','SHH'='red','Group3'='yellow','Group4'='green','n/a (NORM)'='grey')
)

```



### Stage IV: Prepare SJARACNe input (Mouse)
```bash

# ---- 1) Load TF/SIG database (mouse) ----
db.preload(use_level = "gene", use_spe = "mouse", update = FALSE)

# ---- 2) Decide gene ID type & build TF/SIG lists ----
# If fData has gene symbols, use 'external_gene_name'; if it's ENSMUSG..., use 'ensembl_gene_id'.

use_gene_type <- 'external_gene_name'   # 或 'ensembl_gene_id'，根据你的数据来改
use_genes     <- rownames(fData(network.par$net.eset))
use_list      <- get.TF_SIG.list(use_genes, use_gene_type = use_gene_type)

# ---- 3) Choose samples (use all by default) ----
phe         <- pData(network.par$net.eset)
use.samples <- rownames(phe)
stopifnot(length(use.samples) > 0)

# ---- 4) Project name ----
if (is.null(network.par$project.name)) {
  network.par$project.name <- "MyMouseProject"
}
prj.name <- network.par$project.name

# ---- 5) Output directory for SJARACNe inputs ----
network.par$out.dir.SJAR <- file.path(rdir, "SJARACNe")
dir.create(network.par$out.dir.SJAR, showWarnings = FALSE, recursive = TRUE)

# ---- 6) Generate SJARACNe input files ----
SJAracne.prepare(
  eset              = network.par$net.eset,
  use.samples       = use.samples,
  TF_list           = use_list$tf,
  SIG_list          = use_list$sig,
  IQR.thre          = 0.5,
  IQR.loose_thre    = 0.1,
  SJAR.project_name = prj.name,
  SJAR.main_dir     = network.par$out.dir.SJAR
)

# (Optional) Quick peek at created files
cat("[INFO] SJARACNe input dir:\n")
print(file.path(network.par$out.dir.SJAR, prj.name))
print(list.files(file.path(network.par$out.dir.SJAR, prj.name), recursive = TRUE))

```

## Step 5. Run SJARACNe for Network Inference (Server)

We use **SJARACNe** to construct regulatory networks (TF network and signature network).
**Run this on the server after Step 4** has produced: `input.exp`, `tf.txt`, `sig.txt`.\

**Prerequisites (same conda env as** `sjaracne`)
`conda install -c conda-forge numpy pandas networkx cwltool`

**Required inputs**
- `input.exp` → expression matrix (rows = genes, columns = samples; tab-delimited; first column = gene IDs)  
- `tf.txt` → list of transcription factors (one per line, matching IDs in `input.exp`)  
- `sig.txt` → list of signature genes (one per line, matching IDs in `input.exp`)  


```bash
cd /mnt/sda/Public/Project/collabration/AoLab/20250821/3.R_analysis/SJARACNe/MyMouseProject

# 0) Preflight: verify inputs exist
ls -lh input.exp tf.txt sig.txt || { echo "[ERR] Missing input files"; exit 1; }

# 1) Locate the consensus script (used by the CWL workflow)
CONS=$(python - <<'PY'
import os, glob, SJARACNe
base = os.path.dirname(SJARACNe.__file__)
cands = glob.glob(os.path.join(base, "**", "create_consensus_network.py"), recursive=True)
print(cands[0] if cands else "")
PY
)
echo "consensus script: $CONS"

# 2) Clean previous outputs and temp dirs
rm -rf output_tf output_sig /tmp/tmp_tf /tmp/tmp_sig

# 3) Run SJARACNe (bootstrap + consensus)
#    These parameters match the working configuration.
sjaracne local \
  -e "$PWD/input.exp" \
  -g "$PWD/tf.txt" \
  -o output_tf \
  -tmp /tmp/tmp_tf

sjaracne local \
  -e "$PWD/input.exp" \
  -g "$PWD/sig.txt" \
  -o output_sig \
  -tmp /tmp/tmp_sig

# 4) Verify consensus outputs
ls -lh output_tf  | grep -i consensus || echo "[ERR] TF consensus missing"
ls -lh output_sig | grep -i consensus || echo "[ERR] SIG consensus missing"
# Expected: output_*/consensus_network_ncol_.txt

```
**✅ Expected Outputs**
| File / Folder                               | Description                                  |
| ------------------------------------------- | -------------------------------------------- |
| `.../MyMouseProject/output_tf/`             | TF network outputs (bootstrap + consensus)   |
| `.../MyMouseProject/output_tf/consensus_*`  | Consensus TF network edge list (ncol format) |
| `.../MyMouseProject/output_sig/`            | Signature network outputs                    |
| `.../MyMouseProject/output_sig/consensus_*` | Consensus signature network edge list        |



## Step 6. NetBID2 Hidden Driver Estimation
**Tip**: set an environment variable once and forget about editing paths later:
`export PROJECT_ROOT="/mnt/sda/Public/Project/collabration/AoLab/20250821"`
```bash


# ---- Environment ----
suppressPackageStartupMessages({
  library(NetBID2)
  library(Biobase)
})

# ---- Project roots (edit once, or set env PROJECT_ROOT beforehand) ----
project_root <- Sys.getenv("PROJECT_ROOT",
                           "/mnt/sda/Public/Project/collabration/AoLab/20250821")

DIR <- list(
  r       = file.path(project_root, "3.R_analysis"),
  result  = file.path(project_root, "3.R_analysis", "result"),
  sjar    = file.path(project_root, "3.R_analysis", "SJARACNe", "MyMouseProject")
)

project_name <- "KO_vs_WT"

# ---- Working dir ----
dir.create(DIR$r, showWarnings = FALSE, recursive = TRUE)
setwd(DIR$r)

# ---- Locate SJARACNe consensus networks (robust to suffixes) ----
tf_network_file  <- Sys.glob(file.path(DIR$sjar, "output_tf",  "consensus_network_ncol*.txt"))[1]
sig_network_file <- Sys.glob(file.path(DIR$sjar, "output_sig", "consensus_network_ncol*.txt"))[1]

if (is.na(tf_network_file)  || !file.exists(tf_network_file))
  stop("[ERR] TF consensus not found under: ", file.path(DIR$sjar, "output_tf"))
if (is.na(sig_network_file) || !file.exists(sig_network_file))
  stop("[ERR] SIG consensus not found under: ", file.path(DIR$sjar, "output_sig"))

# ---- Optional: backup existing analysis.par ----
if (exists("analysis.par")) {
  save(analysis.par, file = "analysis.par.BAK.RData")
  rm(analysis.par); gc()
}

# ---- Create analysis directory & register networks ----
dir.create(DIR$result, showWarnings = FALSE, recursive = TRUE)
analysis.par <- NetBID.analysis.dir.create(
  project_main_dir = DIR$result,
  project_name     = project_name,
  tf.network.file  = tf_network_file,
  sig.network.file = sig_network_file
)

# ---- Read networks (ncol format from SJARACNe) ----
analysis.par$tf.network  <- get.SJAracne.network(network_file = analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file = analysis.par$sig.network.file)

# ---- Network QC (disable HTML to avoid pandoc requirement) ----
draw.network.QC(analysis.par$tf.network$igraph_obj,
                outdir = analysis.par$out.dir.QC,
                prefix = "TF_net_", html_info_limit = FALSE, generate_html = FALSE)

draw.network.QC(analysis.par$sig.network$igraph_obj,
                outdir = analysis.par$out.dir.QC,
                prefix = "SIG_net_", html_info_limit = TRUE, generate_html = FALSE)

# ---- Merge TF & SIG networks ----
analysis.par$merge.network <- merge_TF_SIG.network(
  TF_network  = analysis.par$tf.network,
  SIG_network = analysis.par$sig.network
)

# ---- Choose ExpressionSet for activity calculation ----
# Prefer the ExpressionSet embedded in the network; otherwise try a saved network.par.
if (!is.null(analysis.par$tf.network$eset) &&
    is(analysis.par$tf.network$eset, "ExpressionSet")) {
  analysis.par$cal.eset <- analysis.par$tf.network$eset
} else {
  cand_rdata <- file.path(DIR$result, "DATA", "network.par.Step.net.RData")
  if (file.exists(cand_rdata)) load(cand_rdata)  # should contain network.par$net.eset
  if (exists("network.par") && is(network.par$net.eset, "ExpressionSet")) {
    analysis.par$cal.eset <- network.par$net.eset
  } else {
    stop("[ERR] No ExpressionSet found for activity calculation: ",
         "provide tf.network$eset or network.par$net.eset.")
  }
}

# ---- Compute driver activity (weighted mean) ----
ac_mat <- cal.Activity(
  target_list = analysis.par$merge.network$target_list,
  cal_mat     = exprs(analysis.par$cal.eset),
  es.method   = "weightedmean"
)

# ---- Wrap activity matrix as ExpressionSet ----
analysis.par$merge.ac.eset <- generate.eset(
  exp_mat         = ac_mat,
  phenotype_info  = pData(analysis.par$cal.eset)[colnames(ac_mat), , drop = FALSE],
  feature_info    = NULL,
  annotation_info = "activity in net-dataset"
)

# ---- QC (optional) ----
draw.eset.QC(
  analysis.par$merge.ac.eset,
  outdir         = analysis.par$out.dir.QC,
  intgroup       = NULL,
  do.logtransform= FALSE,
  prefix         = "AC_",
  pre_define     = c(KO = "red", NC = "blue"),
  emb_plot_type  = "2D.interactive",
  generate_html  = FALSE
)

# ---- Save progress ----
NetBID.saveRData(analysis.par = analysis.par, step = "act-get")

```

<img width="6256" height="4167" alt="image" src="https://github.com/user-attachments/assets/fea832e0-f2b2-42e2-a966-9d0be8072c86" />


## 🧮 Step 7. Differential Expression/Activity (KO vs WT) & Master Table

**Input objects**: uses `analysis.par$cal.eset`（expression）and `analysis.par$merge.ac.eset` (activity) from the previous step.
**Group labels**: edit `grp1_label` / `grp0_label` if your metadata uses other names (e.g., `NC`).

```bash
# ---- Setup containers (no hard-coded paths) ----
analysis.par$DE <- list()
analysis.par$DA <- list()

# ---- Metadata & grouping ----
# 'subgroup_col' must exist in pData(.); it's the column defining your groups.
phe_info <- pData(analysis.par$cal.eset)
subgroup_col <- "group"
if (!subgroup_col %in% colnames(phe_info)) {
  stop(sprintf("Phenotype column not found: %s. Please check pData(analysis.par$cal.eset).", subgroup_col))
}

# Case/control labels (edit if your labels differ; e.g., control = "NC").
grp1_label <- "KO"
grp0_label <- "WT"

# Sample indices per group
G1 <- rownames(phe_info)[phe_info[[subgroup_col]] == grp1_label]
G0 <- rownames(phe_info)[phe_info[[subgroup_col]] == grp0_label]

# Minimal replicate check (recommended ≥ 2 per group)
if (length(G1) < 2 || length(G0) < 2) {
  stop(sprintf("Too few samples per group: %s=%d, %s=%d; need ≥2 per group.",
               grp1_label, length(G1), grp0_label, length(G0)))
}

# Contrast name used in downstream outputs
comp_name <- sprintf("%s.Vs.%s", grp1_label, grp0_label)

# ---- Differential Expression (genes) on expression eSet ----
DE_gene_bid <- getDE.BID.2G(
  eset    = analysis.par$cal.eset,
  G1      = G1, G0 = G0,
  G1_name = grp1_label, G0_name = grp0_label
)

# ---- Differential Activity (drivers) on activity eSet ----
DA_driver_bid <- getDE.BID.2G(
  eset    = analysis.par$merge.ac.eset,
  G1      = G1, G0 = G0,
  G1_name = grp1_label, G0_name = grp0_label
)

######## For a single contrast, do this ###########
## 1) Save results into analysis.par (required for later steps)
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

## 2) Helper: build the structure required by draw.combineDE(), even for a single contrast
make_DElist_for_single <- function(de_df, comp_name,
                                   z_col = "Z-statistics",
                                   display_col = "P.Value") {
  stopifnot(is.data.frame(de_df))
  stopifnot(all(c("ID", z_col, display_col) %in% colnames(de_df)))

  ord  <- order(de_df[[display_col]], na.last = NA)

  # 'combine' = de_df sorted by display_col, with an extra column named by the contrast
  # that stores original rownames (used by draw.combineDE for indexing).
  comb <- de_df[ord, , drop = FALSE]
  comb[[comp_name]] <- rownames(de_df)[ord]
  rownames(comb) <- comb$ID  # ensure rownames align with ID

  out <- list()
  out[[comp_name]] <- de_df   # per-contrast table
  out$combine      <- comb    # combined table
  out
}

## 3) Plot combined DE/DA tables
DE_list <- make_DElist_for_single(analysis.par$DE[[comp_name]], comp_name)
DA_list <- make_DElist_for_single(analysis.par$DA[[comp_name]], comp_name)

pdf_file_DE <- sprintf("%s/%s_DE.pdf", analysis.par$out.dir.PLOT, comp_name)
pdf_file_DA <- sprintf("%s/%s_DA.pdf", analysis.par$out.dir.PLOT, comp_name)

draw.combineDE(DE_list, pdf_file = pdf_file_DE)
draw.combineDE(DA_list, pdf_file = pdf_file_DA)

## ---- V. Generate master table and save ----
# Species DB for annotation; set use_spe = "human" if your data are human.
db.preload(use_level = "gene", use_spe = "mouse", update = FALSE)

all_comp <- names(analysis.par$DE)  # here we have only one contrast

# IMPORTANT: ensure 'tf_sigs' exists; if not, use an empty list:
# tf_sigs <- list()

analysis.par$final_ms_tab <- generate.masterTable(
  use_comp     = all_comp,
  DE           = analysis.par$DE,
  DA           = analysis.par$DA,
  target_list  = analysis.par$merge.network$target_list,
  tf_sigs      = tf_sigs,                 # empty list is OK
  z_col        = "Z-statistics",
  display_col  = c("logFC", "P.Value"),
  main_id_type = "external_gene_name"
)

out_file <- sprintf("%s/%s_ms_tab.xlsx", analysis.par$out.dir.DATA, analysis.par$project.name)

## Optional highlighting (edit markers/colours as needed)
mark_gene <- list(KO = c("GPR160"))
mark_col  <- list(KO = "red", WT = "blue")

out2excel(analysis.par$final_ms_tab, out.xlsx = out_file, mark_gene, mark_col)

## Save RData (essential)
NetBID.saveRData(analysis.par = analysis.par, step = "ms-tab")

cat("\n==== Step 7 completed (KO vs WT) ====\n",
    "Master table: ", out_file, "\n",
    "RData: ", file.path(analysis.par$out.dir.DATA, "analysis.par.Step.ms-tab.RData"), "\n", sep = "")

```


## 🚀 Step 8. Advanced Analysis (Volcano, GSEA, Enrichment, Heatmaps & Network)
This section assumes you have completed **Step 7** and saved **analysis.par** at the ms-tab stage.

**A) Initialize & Folder Layout**
- **Do**: reload objects saved at `ms-tab`, set contrast name, create neat subfolders.

- **Need**: `analysis.par.Step.ms-tab.RData` under your `.../DATA/`.

- **Outputs**: subfolders under `analysis.par$out.dir.PLOT`.
```bash
suppressPackageStartupMessages({
  library(NetBID2); library(Biobase); library(grid)  # grid for gpar()
})

## 1) Contrast name  (edit if needed)
COMP <- "KO.Vs.WT"

## 2) Reload ms-tab step (edit DATA path if needed)
analysis.par <- list()
analysis.par$out.dir.DATA <- "/mnt/sda/Public/Project/collabration/AoLab/20250821/3.R_analysis/result/KO_vs_WT/DATA"
NetBID.loadRData(analysis.par = analysis.par, step = "ms-tab")

## 3) Shorthands
ms_tab   <- analysis.par$final_ms_tab
phe_info <- pData(analysis.par$cal.eset)
exp_mat  <- exprs(analysis.par$cal.eset)
ac_mat   <- exprs(analysis.par$merge.ac.eset)

## 4) Create tidy output buckets
OUTBASE <- analysis.par$out.dir.PLOT
OUT <- list(
  volcano     = file.path(OUTBASE, "01_volcano"),
  heatmap     = file.path(OUTBASE, "02_heatmap"),
  gsea_netbid = file.path(OUTBASE, "03_gsea_netbid"),
  gsea_driver = file.path(OUTBASE, "04_gsea_driver"),
  enrich      = file.path(OUTBASE, "05_enrichment"),
  bubble      = file.path(OUTBASE, "06_bubble"),
  network     = file.path(OUTBASE, "07_network"),
  category    = file.path(OUTBASE, "08_category")
)
invisible(lapply(OUT, dir.create, showWarnings = FALSE, recursive = TRUE))

## helper to compose paths
p <- function(bucket, fname) file.path(OUT[[bucket]], fname)

## Optional: light filter for driver size
ms_tab <- ms_tab[ms_tab$Size >= 30 & ms_tab$Size <= 1000, ]
```


**B) Volcano Plots (DA / DE)**
- **Do**: identify significant drivers/genes and save volcano PDFs.

- **Need**:`ms_tab` from Step A.

- **Outputs**:`01_volcano/volcano_DA_<COMP>.pdf`, `01_volcano/volcano_DE_<COMP>.pdf`.

```bash
sig_driver <- draw.volcanoPlot(
  dat = ms_tab, label_col = "gene_label",
  logFC_col = sprintf("logFC.%s_DA", COMP),
  Pv_col    = sprintf("P.Value.%s_DA", COMP),
  logFC_thre = 0.05, Pv_thre = 0.05, show_label = TRUE, label_cex = 1,
  main = sprintf("Volcano (DA) – %s", COMP),
  pdf_file = p("volcano", sprintf("volcano_DA_%s.pdf", COMP))
)

sig_gene <- draw.volcanoPlot(
  dat = ms_tab, label_col = "geneSymbol",
  logFC_col = sprintf("logFC.%s_DE", COMP),
  Pv_col    = sprintf("P.Value.%s_DE", COMP),
  logFC_thre = 0.06, Pv_thre = 0.08, show_label = TRUE, label_cex = 1,
  main = sprintf("Volcano (DE) – %s", COMP),
  pdf_file = p("volcano", sprintf("volcano_DE_%s.pdf", COMP))
)

driver_list <- rownames(sig_driver)
```

**C) One-shot NetBID GSEA (Drivers → Targets)**
- **Do**: NetBID’s built-in enrichment linking driver activity to target DE ranking.

- **Need**:`analysis.par$DE[[COMP]]`, `driver_list`.

- **Outputs**:`03_gsea_netbid/NetBID_GSEA_*_<COMP>.pdf`.
```bash
DE <- analysis.par$DE[[COMP]]

draw.GSEA.NetBID(
  DE = DE, profile_col = "logFC", profile_trend = "neg2pos", name_col = "ID",
  driver_list = driver_list,
  show_label  = ms_tab[driver_list, "gene_label"],
  driver_DA_Z = ms_tab[driver_list, sprintf("Z.%s_DA", COMP)],
  driver_DE_Z = ms_tab[driver_list, sprintf("Z.%s_DE", COMP)],
  target_list = analysis.par$merge.network$target_list,
  top_driver_number = 30, target_nrow = 2,
  target_col_type = "DE", target_col = "RdBu",
  left_annotation = "high in WT", right_annotation = "high in KO",
  Z_sig_thre = 1.64, profile_sig_thre = 1.64,
  main = COMP,
  pdf_file = p("gsea_netbid", sprintf("NetBID_GSEA_DEprofile_%s.pdf", COMP))
)

draw.GSEA.NetBID(
  DE = DE, profile_col = "t", profile_trend = "pos2neg",
  driver_list = driver_list,
  show_label  = ms_tab[driver_list, "gene_label"],
  driver_DA_Z = ms_tab[driver_list, sprintf("Z.%s_DA", COMP)],
  driver_DE_Z = ms_tab[driver_list, sprintf("Z.%s_DE", COMP)],
  target_list = analysis.par$merge.network$target_list,
  top_driver_number = 30, target_nrow = 2,
  target_col_type = "PN", target_col = "RdBu",
  left_annotation = "high in KO", right_annotation = "high in WT",
  Z_sig_thre = 1.64, profile_sig_thre = 1.64,
  main = COMP,
  pdf_file = p("gsea_netbid", sprintf("NetBID_GSEA_tprofile_%s.pdf", COMP))
)
```
**D) Heatmaps (Expression & Activity)**
- **Do**: visualize expression of top DE genes & activity of top DA drivers.

- **Need**:`sig_gene`, `driver_list`, `exp_mat`, `ac_mat`.

- **Outputs**:`02_heatmap/heatmap_*_<COMP>.pdf`.
```bash
draw.heatmap(
  mat = exp_mat,
  use_genes      = ms_tab[gene_list, "originalID"],
  use_gene_label = ms_tab[gene_list, "geneSymbol"],
  use_samples    = colnames(exp_mat),
  use_sample_label = use_sample_label,
  phenotype_info = phe_info, use_phe = c("group"),
  main = sprintf("Expression (top DE genes) - %s", COMP),  # <-- hyphen
  scale = "row", cluster_rows = TRUE, cluster_columns = TRUE,
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  row_names_gp = gpar(fontsize = 12), pre_define = pre_define,
  pdf_file = p("heatmap", sprintf("heatmap_expression_%s.pdf", COMP))
)

draw.heatmap(
  mat = ac_mat,
  use_genes      = driver_list,
  use_gene_label = ms_tab[driver_list, "gene_label"],
  use_samples    = colnames(ac_mat),
  use_sample_label = use_sample_label,
  phenotype_info = phe_info, use_phe = c("group"),
  main = sprintf("Activity (top DA drivers) - %s", COMP),  # <-- hyphen
  scale = "row", cluster_rows = TRUE, cluster_columns = TRUE,
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  row_names_gp = gpar(fontsize = 12), pre_define = pre_define,
  pdf_file = p("heatmap", sprintf("heatmap_activity_%s.pdf", COMP))
)
```


**E) Functional Enrichment (MSigDB Fisher)**
- **Do**: Hallmark/Reactome/GO/CGP enrichment on up/down DA drivers.

- **Need**:MSigDB mouse gene sets (`gs.preload`),`sig_driver`.

- **Outputs**:`05_enrichment/*.pdf` and `*.xlsx`.
```bash
gs.preload(use_spe = "Mus musculus", update = FALSE)

driver_list_up   <- rownames(sig_driver)[sig_driver[, 2] > 0]
driver_list_down <- rownames(sig_driver)[sig_driver[, 2] < 0]
bg_symbols <- unique(ms_tab[, "geneSymbol"])

res_up <- funcEnrich.Fisher(
  input_list = ms_tab[driver_list_up, "geneSymbol"],
  bg_list    = bg_symbols,
  use_gs     = c("H","CP:REACTOME","BP","CGP"),
  Pv_thre    = 0.1, Pv_adj = "none", min_gs_size = 30, max_gs_size = 500
)
res_down <- funcEnrich.Fisher(
  input_list = ms_tab[driver_list_down, "geneSymbol"],
  bg_list    = bg_symbols,
  use_gs     = c("H","CP:REACTOME","BP","CGP"),
  Pv_thre    = 0.1, Pv_adj = "none", min_gs_size = 30, max_gs_size = 500
)

out2excel(list(UP = res_up, DOWN = res_down),
          out.xlsx = p("enrich", sprintf("funcEnrich_Fisher_%s.xlsx", COMP)))

draw.funcEnrich.bar(res_up,   top_number = 30,
                    main = sprintf("Enrichment (UP) – %s", COMP),
                    pdf_file = p("enrich", sprintf("funcEnrich_bar_UP_%s.pdf", COMP)))
draw.funcEnrich.bar(res_down, top_number = 30,
                    main = sprintf("Enrichment (DOWN) – %s", COMP),
                    pdf_file = p("enrich", sprintf("funcEnrich_bar_DOWN_%s.pdf", COMP)))

draw.funcEnrich.cluster(
  funcEnrich_res = res_up, top_number = 30,
  gs_cex = 1.0, gene_cex = 1.0, pv_cex = 0.9,
  cluster_gs = TRUE, cluster_gene = TRUE, h = 0.95,
  pdf_file = p("enrich", sprintf("funcEnrich_cluster_UP_%s.pdf", COMP))
)
```

**F) Bubble Plots (Targets of DA Drivers)**

- **Do**: bubbleplots of target gene-set enrichments per up/down drivers.

- **Need**: Ensembl↔symbol transfer table.

- **Outputs**: 06_bubble/bubble_*_<COMP>.pdf.

```bash
transfer_tab <- get_IDtransfer2symbol2type(
  from_type = "ensembl_gene_id",
  dataset = "mmusculus_gene_ensembl",  # change to hsapiens_gene_ensembl for human
  use_level = "gene", ignore_version = TRUE
)

draw.bubblePlot(
  driver_list = driver_list_up,
  show_label  = ms_tab[driver_list_up, "gene_label"],
  Z_val       = ms_tab[driver_list_up, sprintf("Z.%s_DA", COMP)],
  driver_type = ms_tab[driver_list_up, "gene_biotype"],
  target_list = analysis.par$merge.network$target_list,
  transfer2symbol2type = transfer_tab,
  bg_list = ms_tab[, "geneSymbol"],
  min_gs_size = 5, max_gs_size = 500,
  use_gs = c("H"), top_geneset_number = 30, top_driver_number = 10,
  main = sprintf("Bubble – UP DA drivers – %s", COMP),
  pdf_file = p("bubble", sprintf("bubble_UP_%s.pdf", COMP))
)

draw.bubblePlot(
  driver_list = driver_list_down,
  show_label  = ms_tab[driver_list_down, "gene_label"],
  Z_val       = ms_tab[driver_list_down, sprintf("Z.%s_DA", COMP)],
  driver_type = ms_tab[driver_list_down, "gene_biotype"],
  target_list = analysis.par$merge.network$target_list,
  transfer2symbol2type = transfer_tab,
  bg_list = ms_tab[, "geneSymbol"],
  min_gs_size = 5, max_gs_size = 500,
  use_gs = c("H"), top_geneset_number = 30, top_driver_number = 10,
  main = sprintf("Bubble – DOWN DA drivers – %s", COMP),
  pdf_file = p("bubble", sprintf("bubble_DOWN_%s.pdf", COMP))
)

```

**G) Deep-Dive a Selected Driver (GSEA + Subnetwork + Category)**

- **Do**: classic GSEA for one driver, starburst subnetwork, and KO/WT category plots.

- **Need**: `driver_list`, `DE`, `target_list`.

- **Outputs**: PDFs in `04_gsea_driver/`,` 07_network/`,` 08_category/`.
```bash
use_driver <- driver_list[1]   # pick the top one (or set manually)

## GSEA (classic)
DE_profile <- DE$`Z-statistics`; names(DE_profile) <- rownames(DE)
use_targets   <- analysis.par$merge.network$target_list[[use_driver]]$target
use_direction <- sign(analysis.par$merge.network$target_list[[use_driver]]$spearman)
annot_txt     <- sprintf("P-value: %s", signif(ms_tab[use_driver, sprintf("P.Value.%s_DA", COMP)], 2))

draw.GSEA(
  rank_profile     = DE_profile,
  use_genes        = use_targets,
  use_direction    = use_direction,
  main             = sprintf("GSEA – %s", ms_tab[use_driver, "gene_label"]),
  left_annotation  = "high in KO", right_annotation = "high in WT",
  annotation       = annot_txt, annotation_cex = 1.2,
  pdf_file         = p("gsea_driver", sprintf("GSEA_one_driver_%s.pdf", COMP))
)

## Subnetwork (single driver)
edge_score <- with(analysis.par$merge.network$target_list[[use_driver]],
                   MI * sign(spearman))
names(edge_score) <- use_targets

draw.targetNet(
  source_label = ms_tab[use_driver, "gene_label"],
  source_z     = ms_tab[use_driver, sprintf("Z.%s_DA", COMP)],
  edge_score   = edge_score,
  label_cex    = 0.4, n_layer = 4, alphabetical_order = TRUE,
  pdf_file     = p("network", sprintf("targetNet_%s.pdf", COMP))
)

## (Optional) Two drivers
if (length(driver_list) >= 2) {
  use_driver2 <- driver_list[2]
  edge_score2 <- with(analysis.par$merge.network$target_list[[use_driver2]],
                      MI * sign(spearman))
  names(edge_score2) <- analysis.par$merge.network$target_list[[use_driver2]]$target
  all_targets <- unique(analysis.par$merge.network$network_dat$target.symbol)

  draw.targetNet.TWO(
    source1_label = ms_tab[use_driver,  "gene_label"], edge_score1 = edge_score,
    source2_label = ms_tab[use_driver2, "gene_label"], edge_score2 = edge_score2,
    source1_z = ms_tab[use_driver,  sprintf("Z.%s_DA", COMP)],
    source2_z = ms_tab[use_driver2, sprintf("Z.%s_DA", COMP)],
    total_possible_target = all_targets, show_test = TRUE, label_cex = 0.2,
    pdf_file = p("network", sprintf("targetNet_two_%s.pdf", COMP))
  )
}

## Category plot (activity & expression by group)
obs_class <- get_obs_label(phe_info, "group")
draw.categoryValue(
  ac_val        = ac_mat[use_driver, ],
  exp_val       = exp_mat[ms_tab[use_driver, "originalID"], ],
  use_obs_class = obs_class,
  class_order   = c("KO", "WT"), class_srt = 30,
  main_ac       = ms_tab[use_driver, "gene_label"],
  main_exp      = ms_tab[use_driver, "geneSymbol"],
  pre_define    = c(KO = "blue", WT = "red"),
  pdf_file      = p("category", sprintf("categoryValue_%s.pdf", COMP))
)
```
