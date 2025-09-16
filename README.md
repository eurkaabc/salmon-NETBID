## 🧬 Salmon-NETBID2 Workflow

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
##

## 📂 Step 0. Organize FASTQ Files
```bash
#
mkdir -p "/mnt/sda/Public/Project/collabration/AoLab/20250821/1.data"
find "/mnt/sda/Public/Project/collabration/AoLab/20250821/rawdata" -type f -name "*.fastq.gz" -exec cp -n {} "/mnt/sda/Public/Project/collabration/AoLab/20250821/1.data/" \;
```

## 🧪 Step 1. Build Salmon Index
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





## 🧪 Step 2. Salmon Quantification (Batch, Paired-End)

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
## ✅ Expected Outputs

After running this step, you should have:

| File / Folder                              | Description                                         |
| ------------------------------------------ | --------------------------------------------------- |
| `2.salmon/<sample>/quant.sf`               | Transcript-level quantification results             |
| `2.salmon/<sample>/aux_info/meta_info.json`| Mapping statistics (used for summary)               |
| `2.salmon/percent_mapped_summary.tsv`      | Summary table (sample name, mapping %, counts, etc.)|



## 📊 Step 3. Prepare Files for R Analysis

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



## 🧩 Build tx2gene.csv mapping file

Depending on whether your quant.sf transcript IDs retain version numbers (e.g., ENSMUST00000193812.1) or not, generate a matching tx2gene.csv.


## A) Keep version numbers (default if `quant.sf` has them):
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


## B) Remove version numbers (if `quant.sf` lacks them):
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


## ✅ Expected Outputs

After this step, you should have:
| File / Folder                | Description                                                                 |
| ---------------------------- | --------------------------------------------------------------------------- |
| `3.R_analysis/quant.sf.zip`  | Archive of all `quant.sf` files (backup/sharing).                           |
| `3.R_analysis/sampleFile`    | List of sample names (one per line).                                        |
| `3.R_analysis/salmon.output` | Mapping of `sample → quant.sf` file path (used by `tximport` / NetBID2).    |
| `3.R_analysis/tx2gene.csv`   | Transcript-to-gene mapping (with or without version numbers, as generated). |




<img width="6256" height="4167" alt="image" src="https://github.com/user-attachments/assets/fea832e0-f2b2-42e2-a966-9d0be8072c86" />





## 

```bash
# ---- Load necessary libraries ----
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(NetBID2)
  library(Biobase)
  library(edgeR)
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

```

## Stage I: Load & Mapping


```bash
# ---- 1) 读入 gene-level eSet ----
eSet_expGene <- load.exp.RNASeq.demoSalmon(
  salmon_dir         = salmon_dir,
  tx2gene            = tx2,            # 或者传文件路径字符串
  use_phenotype_info = meta,
  use_sample_col     = "sample",
  use_design_col     = "group",
  merge_level        = "gene",
  return_type        = "eset"
)
saveRDS(eSet_expGene, file.path(rdir, "Gene.eset.rawcounts.rds"))
cat("[OK] eSet维度: ", paste(dim(eSet_expGene), collapse=" x "), "\n", sep="")

# ---- 2) 基因ID -> 基因名（mgi_symbol 优先）----
transfer_tab <- get_IDtransfer2symbol2type(
  from_type      = "ensembl_gene_id",
  dataset        = "hsapiens_gene_ensembl",   # ← 改成 Human
  use_level      = "gene",
  ignore_version = TRUE
)
to_col <- if ("hgnc_symbol" %in% colnames(transfer_tab)) "hgnc_symbol" else "external_gene_name"

suppressPackageStartupMessages({
  library(Biobase)
  library(NetBID2)
  library(dplyr)
})

## A) 去版本号 + 保证唯一
old_fn   <- featureNames(eSet_expGene)
fn_nov   <- sub("\\.\\d+$","", old_fn)
fn_nov_u <- make.unique(fn_nov, sep = ".dup")
featureNames(eSet_expGene) <- fn_nov_u

## B) 同步 featureData
fd <- data.frame(ensembl_gene_id = fn_nov,
                 row.names = fn_nov_u,
                 stringsAsFactors = FALSE)
featureData(eSet_expGene) <- new("AnnotatedDataFrame", data = fd)

## C) 执行映射与合并
eSet_expGene2 <- update_eset.feature(
  use_eset         = eSet_expGene,
  use_feature_info = transfer_tab,
  from_feature     = "ensembl_gene_id",
  to_feature       = to_col,
  merge_method     = "median"
)

saveRDS(eSet_expGene2, file.path(rdir, "Gene.eset.mapped.rds"))
cat("[OK] 映射后维度: ", paste(dim(eSet_expGene2), collapse=" x "), "\n", sep="")

```


## Stage II: QC + Save
```bash
# 初始化 network.par
network.par <- list()

# 设置 QC 输出目录
network.par$out.dir.QC <- file.path(rdir, "QC")
dir.create(network.par$out.dir.QC, showWarnings = FALSE, recursive = TRUE)

# 把映射好的 eSet 放进去
network.par$net.eset <- eSet_expGene2
all(colnames(exprs(eSet_expGene2)) == meta$sample)

meta <- meta[match(colnames(exprs(eSet_expGene2)), meta$sample), ]

# QC (intgroup 需要对应你的 meta 列，这里假设是 "group")
draw.eset.QC(
  network.par$net.eset,
  outdir          = network.par$out.dir.QC,
  intgroup        = "group",
  do.logtransform = FALSE,
  prefix          = "beforeQC_",
  generate_html   = FALSE
)

# 定义保存目录（比如就放在 rdir 里）
network.par$out.dir.DATA <- rdir
dir.create(network.par$out.dir.DATA, showWarnings = FALSE, recursive = TRUE)

# 保存
NetBID.saveRData(network.par = network.par, step = "exp-load")
```

## Stage III: Normalization######

```bash
# 取表达矩阵
mat <- exprs(network.par$net.eset)

# 去掉低表达基因（90% 样本低于 5%分位数）
choose1 <- apply(mat <= quantile(mat, probs = 0.05), 1, sum) <= ncol(mat) * 0.90
cat("[FILTER] 低表达基因过滤结果:\n")
print(table(choose1))

mat <- mat[choose1, ]

# 更新 eSet
net_eset <- generate.eset(
  exp_mat        = mat,
  phenotype_info = pData(network.par$net.eset)[colnames(mat), ],
  feature_info   = fData(network.par$net.eset)[rownames(mat), ],
  annotation_info= annotation(network.par$net.eset)
)

# 更新 network.par
network.par$net.eset <- net_eset

# QC after normalization
draw.eset.QC(
  network.par$net.eset,
  outdir          = network.par$out.dir.QC,
  intgroup        = "group",   # ⚠️ 这里用你的 meta 中的分组列，比如 "group"
  do.logtransform = FALSE,
  prefix          = "afterQC_",
  generate_html   = FALSE
)

# 保存
NetBID.saveRData(network.par = network.par, step = "exp-QC")


# (Optional) 样本聚类检查


intgroup <- "group"  # 这里同样要对应 meta 的列名
mat <- exprs(network.par$net.eset)

# 用 IQR 过滤高变基因
choose1 <- IQR.filter(exp_mat=mat, use_genes=rownames(mat), thre=0.8)
cat("[FILTER] 高变基因过滤结果:\n")
print(table(choose1))

mat <- mat[choose1, ]

# K-means 聚类 vs 原始标签
pred_label <- draw.emb.kmeans(
  mat       = mat,
  all_k     = NULL,
  obs_label = get_obs_label(pData(network.par$net.eset), intgroup),
  pre_define= c('WNT'='blue','SHH'='red','Group3'='yellow','Group4'='green','n/a (NORM)'='grey')
)
```

```bash
## Stage III: Prepare SJARACNe input (Mouse)############

# 1) 加载数据库 (小鼠 TF/SIG)
db.preload(use_level = 'gene', use_spe = 'mouse', update = FALSE)

# 2) 转换 gene ID → TF/SIG 列表
# ⚠️ 如果 fData 里是基因符号就用 external_gene_name；如果是 ENSMUSG... 就用 ensembl_gene_id
use_gene_type <- 'external_gene_name'   # 或 'ensembl_gene_id'，根据你的数据来改
use_genes     <- rownames(fData(network.par$net.eset))
use_list      <- get.TF_SIG.list(use_genes, use_gene_type = use_gene_type)

# 3) 选择要用的样本
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe)   # 默认用所有样本

# 4) 设置项目名
if (is.null(network.par$project.name)) {
  network.par$project.name <- "MyMouseProject"
}
prj.name <- network.par$project.name

# 5) 设置 SJARACNe 输出目录
network.par$out.dir.SJAR <- file.path(rdir, "SJARACNe")
dir.create(network.par$out.dir.SJAR, showWarnings = FALSE, recursive = TRUE)

# 6) 生成 SJARACNe 输入文件
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
```
## 🔗 Step 4. Run SJARACNe for Network Inference

We use **SJARACNe** to construct regulatory networks (TF network and signature network).  
Make sure you have:

- `input.exp` → expression matrix (rows = genes, columns = samples; tab-delimited; first column = gene IDs)  
- `tf.txt` → list of transcription factors (one per line, matching IDs in `input.exp`)  
- `sig.txt` → list of signature genes (one per line, matching IDs in `input.exp`)  


```bash
cd /mnt/sda/Public/Project/collabration/AoLab/20250821/3.R_analysis/SJARACNe/MyMouseProject
# 1) 确认共识脚本存在且可执行（避免 cwl 最后一步 Permission denied）
CONS=$(python - <<'PY'
import os, glob, SJARACNe
base = os.path.dirname(SJARACNe.__file__)
cands = glob.glob(os.path.join(base, "**", "create_consensus_network.py"), recursive=True)
print(cands[0] if cands else "")
PY
)
echo "consensus script: $CONS"
# 2) 清理旧产物并创建临时目录
rm -rf output_tf /tmp/tmp_tf

# 3) 直接运行 SJARACNe（自动完成 bootstrap + consensus）
#    这些参数就是你跑成功使用的参数
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


# 4) 验证结果是否生成
ls -lh output_tf | grep -i consensus
# 期望看到：output_tf/consensus_network_ncol_.txt
```





