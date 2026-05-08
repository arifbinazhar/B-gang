#!/usr/bin/env bash
# =============================================================================
# RNA-seq Pipeline: De Novo (Trinity) + Reference-Based (Tuxedo Suite)
# =============================================================================
# Tools required:
#   Quality Control : FastQC, Trimmomatic
#   De Novo         : Trinity, RSEM, Bowtie2, edgeR/DESeq2 (via R)
#   Reference-based : HISAT2 (or TopHat2), StringTie, Cufflinks/Cuffdiff
#   General         : SAMtools, featureCounts, R
#
# Directory structure created:
#   $WORKDIR/
#     raw_reads/         - original FASTQ files
#     fastqc_raw/        - QC on raw reads
#     trimmed/           - adapter-trimmed reads
#     fastqc_trimmed/    - QC on trimmed reads
#     denovo/            - Trinity assembly + quantification
#     reference/         - HISAT2 alignment + StringTie quantification
#     logs/              - all log files
# =============================================================================

set -euo pipefail   # exit on error, unset var, pipe failure
IFS=$'\n\t'

# =============================================================================
# 0. USER CONFIGURATION — Edit this section before running
# =============================================================================

WORKDIR="/path/to/project"          # absolute path to your project directory
THREADS=8                           # number of CPU threads

# ---- Reference-based inputs (required for Tuxedo pipeline) ----
GENOME_FA="/path/to/genome.fa"      # reference genome FASTA
GENOME_GTF="/path/to/annotation.gtf" # GTF annotation file
HISAT2_INDEX="$WORKDIR/hisat2_index/genome" # HISAT2 index prefix (built here)

# ---- Adapter trimming (Trimmomatic) ----
ADAPTERS="/path/to/TruSeq3-PE-2.fa" # adapter file bundled with Trimmomatic
LEADING=3
TRAILING=3
MINLEN=36
SLIDINGWINDOW="4:15"

# ---- Sample sheet: space-separated  SAMPLE_ID  READ1  READ2 ----
# For single-end data, set READ2 to "SE"
declare -A SAMPLES
SAMPLES=(
    ["ctrl_rep1"]="/path/to/ctrl_rep1_R1.fastq.gz /path/to/ctrl_rep1_R2.fastq.gz"
    ["ctrl_rep2"]="/path/to/ctrl_rep2_R1.fastq.gz /path/to/ctrl_rep2_R2.fastq.gz"
    ["treat_rep1"]="/path/to/treat_rep1_R1.fastq.gz /path/to/treat_rep1_R2.fastq.gz"
    ["treat_rep2"]="/path/to/treat_rep2_R1.fastq.gz /path/to/treat_rep2_R2.fastq.gz"
)

# ---- Groups for differential expression ----
CTRL_SAMPLES="ctrl_rep1,ctrl_rep2"
TREAT_SAMPLES="treat_rep1,treat_rep2"

# =============================================================================
# 1. SETUP — Create directories and logging
# =============================================================================

mkdir -p "$WORKDIR"/{fastqc_raw,trimmed,fastqc_trimmed,logs}
mkdir -p "$WORKDIR"/denovo/{assembly,quant,diffexp}
mkdir -p "$WORKDIR"/reference/{hisat2_index,alignments,stringtie,counts,diffexp}

LOG_DIR="$WORKDIR/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
MASTER_LOG="$LOG_DIR/pipeline_${TIMESTAMP}.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$MASTER_LOG"
}

check_tool() {
    command -v "$1" &>/dev/null || { log "ERROR: $1 not found in PATH. Please install it."; exit 1; }
}

log "========== RNA-seq Pipeline Started =========="
log "Working directory : $WORKDIR"
log "Threads           : $THREADS"

# --- Verify all required tools are installed ---
log "Checking required tools..."
for tool in fastqc trimmomatic Trinity bowtie2 rsem-calculate-expression \
            hisat2 samtools stringtie featureCounts Rscript; do
    check_tool "$tool"
done
log "All tools found."

# =============================================================================
# 2. QUALITY CONTROL ON RAW READS
# =============================================================================

log "---------- Step 2: FastQC on raw reads ----------"

for SAMPLE in "${!SAMPLES[@]}"; do
    read -r R1 R2 <<< "${SAMPLES[$SAMPLE]}"
    if [[ "$R2" == "SE" ]]; then
        fastqc -t "$THREADS" -o "$WORKDIR/fastqc_raw" "$R1" \
            &>> "$LOG_DIR/${SAMPLE}_fastqc_raw.log"
    else
        fastqc -t "$THREADS" -o "$WORKDIR/fastqc_raw" "$R1" "$R2" \
            &>> "$LOG_DIR/${SAMPLE}_fastqc_raw.log"
    fi
    log "  FastQC raw done: $SAMPLE"
done

# =============================================================================
# 3. ADAPTER TRIMMING WITH TRIMMOMATIC
# =============================================================================

log "---------- Step 3: Trimmomatic adapter trimming ----------"

for SAMPLE in "${!SAMPLES[@]}"; do
    read -r R1 R2 <<< "${SAMPLES[$SAMPLE]}"
    TRIM_DIR="$WORKDIR/trimmed"

    if [[ "$R2" == "SE" ]]; then
        # Single-end trimming
        trimmomatic SE -threads "$THREADS" \
            "$R1" \
            "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10 \
            LEADING:"$LEADING" TRAILING:"$TRAILING" \
            SLIDINGWINDOW:"$SLIDINGWINDOW" MINLEN:"$MINLEN" \
            &>> "$LOG_DIR/${SAMPLE}_trimmomatic.log"
    else
        # Paired-end trimming
        trimmomatic PE -threads "$THREADS" \
            "$R1" "$R2" \
            "$TRIM_DIR/${SAMPLE}_R1_paired.fastq.gz" \
            "$TRIM_DIR/${SAMPLE}_R1_unpaired.fastq.gz" \
            "$TRIM_DIR/${SAMPLE}_R2_paired.fastq.gz" \
            "$TRIM_DIR/${SAMPLE}_R2_unpaired.fastq.gz" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10 \
            LEADING:"$LEADING" TRAILING:"$TRAILING" \
            SLIDINGWINDOW:"$SLIDINGWINDOW" MINLEN:"$MINLEN" \
            &>> "$LOG_DIR/${SAMPLE}_trimmomatic.log"
    fi
    log "  Trimmomatic done: $SAMPLE"
done

# ---- FastQC on trimmed reads ----
log "  Running FastQC on trimmed reads..."
for SAMPLE in "${!SAMPLES[@]}"; do
    read -r R1 R2 <<< "${SAMPLES[$SAMPLE]}"
    TRIM_DIR="$WORKDIR/trimmed"

    if [[ "$R2" == "SE" ]]; then
        fastqc -t "$THREADS" -o "$WORKDIR/fastqc_trimmed" \
            "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz" \
            &>> "$LOG_DIR/${SAMPLE}_fastqc_trimmed.log"
    else
        fastqc -t "$THREADS" -o "$WORKDIR/fastqc_trimmed" \
            "$TRIM_DIR/${SAMPLE}_R1_paired.fastq.gz" \
            "$TRIM_DIR/${SAMPLE}_R2_paired.fastq.gz" \
            &>> "$LOG_DIR/${SAMPLE}_fastqc_trimmed.log"
    fi
done
log "FastQC on trimmed reads complete."

# =============================================================================
# ===========================================================================
#                PIPELINE A: DE NOVO ASSEMBLY (TRINITY + RSEM)
# ===========================================================================
# =============================================================================

log "========== PIPELINE A: De Novo (Trinity) =========="

# ---------------------------------------------------------------------------
# A1. Collect trimmed reads into comma-separated lists for Trinity
# ---------------------------------------------------------------------------

LEFT_READS=""
RIGHT_READS=""
SINGLE_READS=""
IS_PAIRED=true

for SAMPLE in "${!SAMPLES[@]}"; do
    read -r R1 R2 <<< "${SAMPLES[$SAMPLE]}"
    TRIM_DIR="$WORKDIR/trimmed"

    if [[ "$R2" == "SE" ]]; then
        IS_PAIRED=false
        SINGLE_READS+="$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz,"
    else
        LEFT_READS+="$TRIM_DIR/${SAMPLE}_R1_paired.fastq.gz,"
        RIGHT_READS+="$TRIM_DIR/${SAMPLE}_R2_paired.fastq.gz,"
    fi
done

# Strip trailing commas
LEFT_READS="${LEFT_READS%,}"
RIGHT_READS="${RIGHT_READS%,}"
SINGLE_READS="${SINGLE_READS%,}"

# ---------------------------------------------------------------------------
# A2. Trinity de novo transcriptome assembly
# ---------------------------------------------------------------------------

log "---------- Step A2: Trinity Assembly ----------"

TRINITY_OUT="$WORKDIR/denovo/assembly"

if [[ "$IS_PAIRED" == true ]]; then
    Trinity \
        --seqType fq \
        --max_memory 50G \
        --left  "$LEFT_READS" \
        --right "$RIGHT_READS" \
        --CPU "$THREADS" \
        --output "$TRINITY_OUT" \
        &>> "$LOG_DIR/trinity_assembly.log"
else
    Trinity \
        --seqType fq \
        --max_memory 50G \
        --single "$SINGLE_READS" \
        --CPU "$THREADS" \
        --output "$TRINITY_OUT" \
        &>> "$LOG_DIR/trinity_assembly.log"
fi

TRINITY_FASTA="$TRINITY_OUT/Trinity.fasta"
log "Trinity assembly complete. Transcripts: $TRINITY_FASTA"

# --- Basic assembly stats ---
$TRINITY_HOME/util/TrinityStats.pl "$TRINITY_FASTA" \
    > "$LOG_DIR/trinity_assembly_stats.txt" 2>&1
log "Assembly statistics written to $LOG_DIR/trinity_assembly_stats.txt"

# ---------------------------------------------------------------------------
# A3. Build RSEM reference from Trinity assembly
# ---------------------------------------------------------------------------

log "---------- Step A3: RSEM Prepare Reference ----------"

RSEM_REF="$WORKDIR/denovo/quant/rsem_ref/trinity_rsem"
mkdir -p "$(dirname $RSEM_REF)"

$TRINITY_HOME/util/align_and_estimate_abundance.pl \
    --transcripts "$TRINITY_FASTA" \
    --est_method RSEM \
    --aln_method bowtie2 \
    --trinity_mode \
    --prep_reference \
    --output_dir "$(dirname $RSEM_REF)" \
    &>> "$LOG_DIR/rsem_prep_reference.log"

log "RSEM reference built."

# ---------------------------------------------------------------------------
# A4. Per-sample abundance estimation with RSEM
# ---------------------------------------------------------------------------

log "---------- Step A4: RSEM Abundance Estimation ----------"

RSEM_QUANT_DIR="$WORKDIR/denovo/quant"

for SAMPLE in "${!SAMPLES[@]}"; do
    read -r R1 R2 <<< "${SAMPLES[$SAMPLE]}"
    TRIM_DIR="$WORKDIR/trimmed"
    OUT_DIR="$RSEM_QUANT_DIR/$SAMPLE"
    mkdir -p "$OUT_DIR"

    if [[ "$R2" == "SE" ]]; then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl \
            --transcripts "$TRINITY_FASTA" \
            --seqType fq \
            --single "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz" \
            --est_method RSEM \
            --aln_method bowtie2 \
            --trinity_mode \
            --output_dir "$OUT_DIR" \
            --thread_count "$THREADS" \
            &>> "$LOG_DIR/${SAMPLE}_rsem.log"
    else
        $TRINITY_HOME/util/align_and_estimate_abundance.pl \
            --transcripts "$TRINITY_FASTA" \
            --seqType fq \
            --left  "$TRIM_DIR/${SAMPLE}_R1_paired.fastq.gz" \
            --right "$TRIM_DIR/${SAMPLE}_R2_paired.fastq.gz" \
            --est_method RSEM \
            --aln_method bowtie2 \
            --trinity_mode \
            --output_dir "$OUT_DIR" \
            --thread_count "$THREADS" \
            &>> "$LOG_DIR/${SAMPLE}_rsem.log"
    fi
    log "  RSEM done: $SAMPLE"
done

# ---------------------------------------------------------------------------
# A5. Build transcript expression matrix
# ---------------------------------------------------------------------------

log "---------- Step A5: Build Expression Matrix ----------"

# Gather isoform.results files
ISOFORM_FILES=""
for SAMPLE in "${!SAMPLES[@]}"; do
    ISOFORM_FILES+="$RSEM_QUANT_DIR/$SAMPLE/RSEM.isoforms.results "
done

$TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
    --est_method RSEM \
    --gene_trans_map "$TRINITY_OUT/Trinity.fasta.gene_trans_map" \
    --name_sample_by_basedir \
    --out_prefix "$WORKDIR/denovo/quant/trinity_matrix" \
    $ISOFORM_FILES \
    &>> "$LOG_DIR/expression_matrix.log"

log "Expression matrix built at $WORKDIR/denovo/quant/trinity_matrix.*"

# ---------------------------------------------------------------------------
# A6. Differential Expression with edgeR (via Trinity helper script)
# ---------------------------------------------------------------------------

log "---------- Step A6: Differential Expression (edgeR) ----------"

DE_DIR="$WORKDIR/denovo/diffexp"

$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
    --matrix "$WORKDIR/denovo/quant/trinity_matrix.gene.counts.matrix" \
    --method edgeR \
    --samples_file <(printf "ctrl\t%s\n" ${CTRL_SAMPLES//,/$'\n'ctrl\t}; \
                    printf "treat\t%s\n" ${TREAT_SAMPLES//,/$'\n'treat\t}) \
    --output "$DE_DIR" \
    &>> "$LOG_DIR/edgeR_denovo.log"

log "De novo DE analysis complete. Results in $DE_DIR"

# ---------------------------------------------------------------------------
# A7. Extract significantly DE genes and generate heatmap
# ---------------------------------------------------------------------------

log "---------- Step A7: DE Gene Extraction & Heatmap ----------"

cd "$DE_DIR"
$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl \
    --matrix "$WORKDIR/denovo/quant/trinity_matrix.gene.TMM.EXPR.matrix" \
    -P 0.05 \
    -C 2 \
    &>> "$LOG_DIR/de_analysis_denovo.log"
cd - &>/dev/null

log "De Novo Pipeline (Trinity + RSEM + edgeR) COMPLETE."

# =============================================================================
# ===========================================================================
#         PIPELINE B: REFERENCE-BASED RNA-seq (TUXEDO / HISAT2+StringTie)
# ===========================================================================
# =============================================================================

log "========== PIPELINE B: Reference-Based (Tuxedo/HISAT2+StringTie) =========="

# ---------------------------------------------------------------------------
# B1. Build HISAT2 genome index (skip if already built)
# ---------------------------------------------------------------------------

log "---------- Step B1: Build HISAT2 Index ----------"

HISAT2_IDX_DIR="$WORKDIR/reference/hisat2_index"
mkdir -p "$HISAT2_IDX_DIR"

if [[ ! -f "${HISAT2_IDX_DIR}/genome.1.ht2" ]]; then
    # Extract splice sites and exons from GTF for better alignment
    hisat2_extract_splice_sites.py "$GENOME_GTF" \
        > "$HISAT2_IDX_DIR/splicesites.txt" 2>> "$LOG_DIR/hisat2_index.log"

    hisat2_extract_exons.py "$GENOME_GTF" \
        > "$HISAT2_IDX_DIR/exons.txt" 2>> "$LOG_DIR/hisat2_index.log"

    hisat2-build \
        -p "$THREADS" \
        --ss "$HISAT2_IDX_DIR/splicesites.txt" \
        --exon "$HISAT2_IDX_DIR/exons.txt" \
        "$GENOME_FA" \
        "$HISAT2_IDX_DIR/genome" \
        &>> "$LOG_DIR/hisat2_index.log"
    log "HISAT2 index built."
else
    log "HISAT2 index already exists, skipping."
fi

HISAT2_INDEX="$HISAT2_IDX_DIR/genome"

# ---------------------------------------------------------------------------
# B2. Align each sample to reference with HISAT2
# ---------------------------------------------------------------------------

log "---------- Step B2: HISAT2 Alignment ----------"

ALN_DIR="$WORKDIR/reference/alignments"

for SAMPLE in "${!SAMPLES[@]}"; do
    read -r R1 R2 <<< "${SAMPLES[$SAMPLE]}"
    TRIM_DIR="$WORKDIR/trimmed"
    SAMPLE_BAM="$ALN_DIR/${SAMPLE}.sorted.bam"

    log "  Aligning: $SAMPLE"

    if [[ "$R2" == "SE" ]]; then
        hisat2 \
            -p "$THREADS" \
            --dta \
            -x "$HISAT2_INDEX" \
            -U "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz" \
            --rg-id "$SAMPLE" \
            --rg "SM:$SAMPLE" \
            --summary-file "$LOG_DIR/${SAMPLE}_hisat2_summary.txt" \
            2>> "$LOG_DIR/${SAMPLE}_hisat2.log" \
        | samtools sort -@ "$THREADS" -o "$SAMPLE_BAM"
    else
        hisat2 \
            -p "$THREADS" \
            --dta \
            -x "$HISAT2_INDEX" \
            -1 "$TRIM_DIR/${SAMPLE}_R1_paired.fastq.gz" \
            -2 "$TRIM_DIR/${SAMPLE}_R2_paired.fastq.gz" \
            --rg-id "$SAMPLE" \
            --rg "SM:$SAMPLE" \
            --summary-file "$LOG_DIR/${SAMPLE}_hisat2_summary.txt" \
            2>> "$LOG_DIR/${SAMPLE}_hisat2.log" \
        | samtools sort -@ "$THREADS" -o "$SAMPLE_BAM"
    fi

    samtools index "$SAMPLE_BAM" &>> "$LOG_DIR/${SAMPLE}_samtools.log"

    # Alignment stats
    samtools flagstat "$SAMPLE_BAM" > "$LOG_DIR/${SAMPLE}_flagstat.txt" 2>&1
    log "  Alignment done: $SAMPLE — $(grep 'mapped (' $LOG_DIR/${SAMPLE}_flagstat.txt | head -1)"
done

# ---------------------------------------------------------------------------
# B3. Transcript assembly with StringTie (per sample)
# ---------------------------------------------------------------------------

log "---------- Step B3: StringTie Per-Sample Assembly ----------"

STRINGTIE_DIR="$WORKDIR/reference/stringtie"
mkdir -p "$STRINGTIE_DIR"

for SAMPLE in "${!SAMPLES[@]}"; do
    SAMPLE_BAM="$ALN_DIR/${SAMPLE}.sorted.bam"
    OUT_GTF="$STRINGTIE_DIR/${SAMPLE}.gtf"

    stringtie "$SAMPLE_BAM" \
        -G "$GENOME_GTF" \
        -o "$OUT_GTF" \
        -p "$THREADS" \
        -l "$SAMPLE" \
        &>> "$LOG_DIR/${SAMPLE}_stringtie_assembly.log"

    log "  StringTie assembled: $SAMPLE"
done

# ---------------------------------------------------------------------------
# B4. StringTie merge — produce unified transcriptome
# ---------------------------------------------------------------------------

log "---------- Step B4: StringTie Merge ----------"

# Create a list of per-sample GTFs
find "$STRINGTIE_DIR" -name "*.gtf" \
    ! -name "merged*" \
    ! -name "*.ballgown*" \
    > "$STRINGTIE_DIR/gtf_list.txt"

MERGED_GTF="$STRINGTIE_DIR/merged.gtf"

stringtie --merge \
    -G "$GENOME_GTF" \
    -o "$MERGED_GTF" \
    "$STRINGTIE_DIR/gtf_list.txt" \
    &>> "$LOG_DIR/stringtie_merge.log"

log "Merged GTF: $MERGED_GTF"

# Assess how well merged transcriptome matches known annotation
gffcompare \
    -r "$GENOME_GTF" \
    -G \
    -o "$STRINGTIE_DIR/gffcmp" \
    "$MERGED_GTF" \
    &>> "$LOG_DIR/gffcompare.log"

log "gffcompare annotation comparison done."

# ---------------------------------------------------------------------------
# B5. StringTie re-quantification with merged GTF (Ballgown mode)
# ---------------------------------------------------------------------------

log "---------- Step B5: StringTie Re-quantification ----------"

BALLGOWN_DIR="$WORKDIR/reference/stringtie/ballgown"
mkdir -p "$BALLGOWN_DIR"

for SAMPLE in "${!SAMPLES[@]}"; do
    SAMPLE_BAM="$ALN_DIR/${SAMPLE}.sorted.bam"
    OUT_DIR="$BALLGOWN_DIR/$SAMPLE"
    mkdir -p "$OUT_DIR"

    stringtie "$SAMPLE_BAM" \
        -G "$MERGED_GTF" \
        -o "$OUT_DIR/${SAMPLE}.gtf" \
        -p "$THREADS" \
        -e \
        -B \
        &>> "$LOG_DIR/${SAMPLE}_stringtie_quant.log"

    log "  StringTie re-quant done: $SAMPLE"
done

# ---------------------------------------------------------------------------
# B6. featureCounts — raw count matrix (for DESeq2)
# ---------------------------------------------------------------------------

log "---------- Step B6: featureCounts Raw Count Matrix ----------"

COUNTS_DIR="$WORKDIR/reference/counts"
BAM_FILES=""
for SAMPLE in "${!SAMPLES[@]}"; do
    BAM_FILES+="$ALN_DIR/${SAMPLE}.sorted.bam "
done

featureCounts \
    -T "$THREADS" \
    -a "$GENOME_GTF" \
    -o "$COUNTS_DIR/gene_counts.txt" \
    -p \
    --countReadPairs \
    -B \
    -C \
    $BAM_FILES \
    &>> "$LOG_DIR/featureCounts.log"

log "featureCounts matrix: $COUNTS_DIR/gene_counts.txt"

# ---------------------------------------------------------------------------
# B7. Differential Expression — DESeq2 via R
# ---------------------------------------------------------------------------

log "---------- Step B7: DESeq2 Differential Expression ----------"

DE_DIR_REF="$WORKDIR/reference/diffexp"

# Build sample metadata table
META_FILE="$DE_DIR_REF/sample_metadata.csv"
{
    echo "sample,condition"
    for S in ${CTRL_SAMPLES//,/ };   do echo "$S,ctrl";  done
    for S in ${TREAT_SAMPLES//,/ };  do echo "$S,treat"; done
} > "$META_FILE"

# Write DESeq2 R script inline
DESEQ2_SCRIPT="$DE_DIR_REF/deseq2_analysis.R"

cat > "$DESEQ2_SCRIPT" << 'RSCRIPT'
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

args    <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
meta_file   <- args[2]
out_dir     <- args[3]

# --- Load count matrix ---
counts_raw <- read.table(counts_file, header = TRUE, row.names = 1,
                         sep = "\t", comment.char = "#")

# featureCounts: columns 1-5 are metadata; BAM counts start at col 6
counts_mat <- as.matrix(counts_raw[, 6:ncol(counts_raw)])
colnames(counts_mat) <- sub(".sorted.bam", "", basename(colnames(counts_mat)))

# --- Sample metadata ---
meta <- read.csv(meta_file, row.names = 1)
meta$condition <- factor(meta$condition, levels = c("ctrl", "treat"))
meta <- meta[colnames(counts_mat), , drop = FALSE]

# --- DESeq2 ---
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = meta,
                              design    = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]   # pre-filter low counts
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "treat", "ctrl"),
               alpha = 0.05)
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), file = file.path(out_dir, "DESeq2_results.csv"))

# --- Significant genes ---
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 1)
write.csv(as.data.frame(sig), file = file.path(out_dir, "DESeq2_significant.csv"))
message(nrow(sig), " significant DE genes (padj<0.05, |LFC|>=1)")

# --- PCA plot ---
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
png(file.path(out_dir, "PCA_plot.png"), width = 800, height = 600)
print(
    ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
        geom_point(size = 4) +
        geom_text(vjust = -0.8, size = 3) +
        xlab(paste0("PC1: ", pct_var[1], "% variance")) +
        ylab(paste0("PC2: ", pct_var[2], "% variance")) +
        theme_classic(base_size = 14) +
        ggtitle("PCA of samples (VST-normalized)")
)
dev.off()

# --- Heatmap of top 50 DE genes ---
top50 <- head(rownames(sig), 50)
if (length(top50) >= 2) {
    mat <- assay(vsd)[top50, ]
    mat <- mat - rowMeans(mat)
    anno <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])
    png(file.path(out_dir, "heatmap_top50.png"), width = 900, height = 1000)
    pheatmap(mat,
             annotation_col   = anno,
             show_rownames    = TRUE,
             show_colnames    = TRUE,
             cluster_rows     = TRUE,
             cluster_cols     = TRUE,
             fontsize_row     = 8,
             main             = "Top 50 DE Genes (VST, row-centered)")
    dev.off()
}

# --- Volcano plot ---
res_df <- as.data.frame(res)
res_df$sig <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 &
                     abs(res_df$log2FoldChange) >= 1, "Significant", "NS")
png(file.path(out_dir, "volcano_plot.png"), width = 900, height = 700)
print(
    ggplot(res_df, aes(log2FoldChange, -log10(padj), color = sig)) +
        geom_point(alpha = 0.5, size = 1.5) +
        scale_color_manual(values = c("Significant" = "red", "NS" = "grey60")) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        theme_classic(base_size = 14) +
        labs(title = "Volcano Plot (treat vs ctrl)",
             x = "log2 Fold Change", y = "-log10(adjusted p-value)")
)
dev.off()

message("DESeq2 analysis complete.")
RSCRIPT

Rscript "$DESEQ2_SCRIPT" \
    "$COUNTS_DIR/gene_counts.txt" \
    "$META_FILE" \
    "$DE_DIR_REF" \
    &>> "$LOG_DIR/deseq2.log"

log "DESeq2 analysis complete. Results in $DE_DIR_REF"

# ---------------------------------------------------------------------------
# B8. (Optional) Ballgown — transcript-level DE using StringTie output
# ---------------------------------------------------------------------------

log "---------- Step B8: Ballgown Transcript-level DE ----------"

BALLGOWN_SCRIPT="$DE_DIR_REF/ballgown_analysis.R"

cat > "$BALLGOWN_SCRIPT" << 'RSCRIPT'
suppressPackageStartupMessages({
    library(ballgown)
    library(RSkittleBrewer)
    library(genefilter)
    library(dplyr)
})

args         <- commandArgs(trailingOnly = TRUE)
bg_dir       <- args[1]   # path to ballgown per-sample dirs
meta_file    <- args[2]
out_dir      <- args[3]

meta <- read.csv(meta_file)

# Load ballgown object
bg <- ballgown(dataDir = bg_dir, samplePattern = "", pData = meta,
               verbose = FALSE)

# Filter low-variance transcripts
bg_filt <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset = TRUE)

# Transcript-level DE
results_tx <- stattest(bg_filt, feature = "transcript",
                       covariate = "condition",
                       adjustvars = NULL, getFC = TRUE, meas = "FPKM")
results_tx <- results_tx %>% arrange(qval)
write.csv(results_tx, file = file.path(out_dir, "ballgown_transcript_DE.csv"),
          row.names = FALSE)

# Gene-level DE
results_gene <- stattest(bg_filt, feature = "gene",
                         covariate = "condition",
                         adjustvars = NULL, getFC = TRUE, meas = "FPKM")
results_gene <- results_gene %>% arrange(qval)
write.csv(results_gene, file = file.path(out_dir, "ballgown_gene_DE.csv"),
          row.names = FALSE)

sig_tx   <- subset(results_tx,   results_tx$qval   < 0.05)
sig_gene <- subset(results_gene, results_gene$qval < 0.05)

message("Significant transcripts: ", nrow(sig_tx))
message("Significant genes      : ", nrow(sig_gene))
message("Ballgown analysis complete.")
RSCRIPT

Rscript "$BALLGOWN_SCRIPT" \
    "$BALLGOWN_DIR" \
    "$META_FILE" \
    "$DE_DIR_REF" \
    &>> "$LOG_DIR/ballgown.log"

log "Ballgown analysis complete."

# =============================================================================
# 4. FINAL SUMMARY
# =============================================================================

log "========== Pipeline Summary =========="
log ""
log "── De Novo (Trinity + RSEM + edgeR) ──────────────────────────────────"
log "  Assembly    : $TRINITY_FASTA"
log "  Quant dir   : $WORKDIR/denovo/quant"
log "  DE results  : $WORKDIR/denovo/diffexp"
log ""
log "── Reference-Based (HISAT2 + StringTie + DESeq2 + Ballgown) ──────────"
log "  HISAT2 BAMs : $ALN_DIR"
log "  Merged GTF  : $MERGED_GTF"
log "  Counts      : $COUNTS_DIR/gene_counts.txt"
log "  DESeq2 out  : $DE_DIR_REF/DESeq2_results.csv"
log "  Ballgown    : $DE_DIR_REF/ballgown_gene_DE.csv"
log "  Plots       : PCA, Heatmap, Volcano in $DE_DIR_REF"
log ""
log "  All logs    : $LOG_DIR"
log "  Master log  : $MASTER_LOG"
log ""
log "========== Pipeline Finished Successfully =========="