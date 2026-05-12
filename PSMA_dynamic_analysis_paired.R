#!/usr/bin/env Rscript

# =============================================================================
# PSMA DYNAMICS ANALYSIS — PAIRED PATIENTS (INITIAL → RECURRENT)
# Auto-detects patients with both timepoints, then runs the full analysis.
#
# Expected directory layout (outputs are created automatically):
#
#   <project_root>/                          ← where this script lives
#     PSMA_dynamic_analysis_paired.R
#       TCGA_GBM_IDHwt_with_RNAseq/
#         rna_seq_raw_counts_initial_recurrent_filtered.csv
#         TCGA_GBM_clinical_initial_recurrent.csv
#     PSMA_Paired_Dynamics/                ← created automatically
#         01_PSMA_Evolution_Overview.png
#         02_Volcano_DGE.png
#         03_Heatmap_Top50_DEGs.png
#         04_Biological_Scores_TTest.png
#         05_KM_Survival_Analysis.png
#         06_GO_Enrichment_Increased_Group.csv
#         06_GO_Dotplot.png
#         07_GO_Barplot.png
#         DGE_Full_Results.csv
#         paired_patients_report.txt
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Portable path detection
# -----------------------------------------------------------------------------
# Works in RStudio (interactive / sourced), Rscript --file=, and as a fallback
# to the current working directory.
get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0)
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}

script_dir    <- get_script_dir()
message("Script directory: ", script_dir)

data_dir      <- file.path(script_dir, "TCGA_GBM_IDHwt_with_RNAseq")
counts_file   <- file.path(data_dir, "rna_seq_raw_counts_initial_recurrent_filtered.csv")
clinical_file <- file.path(data_dir, "TCGA_GBM_clinical_initial_recurrent.csv")
output_dir    <- file.path(script_dir, "PSMA_Paired_Dynamics")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

folh1_id <- "ENSG00000086205"   # Ensembl ID for FOLH1 / PSMA

# -----------------------------------------------------------------------------
# 1. Package installation and loading
# -----------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "enrichplot")
cran_pkgs <- c("dplyr", "tibble", "ggplot2", "ggrepel",
               "survival", "survminer", "msigdbr",
               "patchwork", "pheatmap", "ggpubr")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(survival)
  library(survminer)
  library(msigdbr)
  library(patchwork)
  library(pheatmap)
  library(ggpubr)       # for stat_compare_means()
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

# -----------------------------------------------------------------------------
# 2. Data loading
# -----------------------------------------------------------------------------
message("Step 1: Loading raw counts...")

if (!file.exists(counts_file))
  stop("Counts file not found. Check that 'counts_file' points to the right location:\n  ",
       counts_file)

raw_df     <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts_raw <- if (any(grepl("TCGA", colnames(raw_df)))) as.matrix(raw_df) else t(as.matrix(raw_df))

# Standardise TCGA barcodes ("." → "-") and strip Ensembl version suffixes
colnames(counts_raw) <- gsub("\\.", "-", colnames(counts_raw))
rownames(counts_raw) <- gsub("\\.\\d+$", "", rownames(counts_raw))

# -----------------------------------------------------------------------------
# 3. Auto-detect paired patients (replaces the hardcoded patient list)
# -----------------------------------------------------------------------------
message("Step 2: Identifying patients with BOTH initial and recurrent samples...")

all_barcodes <- colnames(counts_raw)

df_barcodes <- data.frame(
  sample_barcode = all_barcodes,
  patient_id     = substr(all_barcodes, 1, 12),
  type_code      = substr(all_barcodes, 14, 15),
  stringsAsFactors = FALSE
)

initial_pts   <- df_barcodes %>% filter(type_code == "01") %>% pull(patient_id) %>% unique()
recurrent_pts <- df_barcodes %>% filter(type_code == "02") %>% pull(patient_id) %>% unique()

paired_pts    <- intersect(initial_pts, recurrent_pts)
only_initial  <- setdiff(initial_pts,   recurrent_pts)
only_recurrent <- setdiff(recurrent_pts, initial_pts)

# Print overlap report to console
cat("====================================================\n")
cat("         PATIENT OVERLAP REPORT (TCGA)\n")
cat("====================================================\n")
cat("Total barcodes read      :", nrow(df_barcodes),   "\n\n")
cat("Initial patients   (01)  :", length(initial_pts),   "\n")
cat("Recurrent patients (02)  :", length(recurrent_pts), "\n")
cat("----------------------------------------------------\n")
cat("Paired (both timepoints) :", length(paired_pts),    "\n")
cat("Initial only             :", length(only_initial),  "\n")
cat("Recurrent only           :", length(only_recurrent),"\n")
cat("====================================================\n")
if (length(paired_pts) > 0) {
  cat("\nPaired patient IDs:\n")
  print(paired_pts)
}

# Save report to file
writeLines(c(
  "PATIENT OVERLAP REPORT",
  paste("Total barcodes       :", nrow(df_barcodes)),
  paste("Initial (01)         :", length(initial_pts)),
  paste("Recurrent (02)       :", length(recurrent_pts)),
  paste("Paired               :", length(paired_pts)),
  paste("Initial only         :", length(only_initial)),
  paste("Recurrent only       :", length(only_recurrent)),
  "",
  "Paired patient IDs:",
  paired_pts
), file.path(output_dir, "paired_patients_report.txt"))

if (length(paired_pts) < 2)
  stop("At least 2 paired patients are required to run the dynamic analysis. ",
       "Only ", length(paired_pts), " found.")

message("  -> ", length(paired_pts), " paired patients identified. Proceeding with analysis.")

# -----------------------------------------------------------------------------
# 4. Build paired metadata and count matrix
# -----------------------------------------------------------------------------
message("Step 3: Building paired metadata and count matrix...")

meta_paired <- data.frame(sample_barcode = colnames(counts_raw)) %>%
  mutate(
    patient_id = substr(sample_barcode, 1, 12),
    type       = ifelse(substr(sample_barcode, 14, 15) == "01", "Initial", "Recurrent")
  ) %>%
  filter(patient_id %in% paired_pts) %>%
  distinct(patient_id, type, .keep_all = TRUE) %>%   # one sample per timepoint per patient
  arrange(patient_id, type)

counts_paired <- counts_raw[, meta_paired$sample_barcode]
storage.mode(counts_paired) <- "integer"
counts_paired[is.na(counts_paired)] <- 0L
counts_paired <- counts_paired[rowSums(counts_paired) > 0, ]

if (!folh1_id %in% rownames(counts_paired))
  stop("FOLH1 (", folh1_id, ") is absent from the count matrix after filtering. ",
       "Check that the correct gene ID is used.")

# -----------------------------------------------------------------------------
# 5. DESeq2 normalisation and PSMA trend classification
# -----------------------------------------------------------------------------
message("Step 4: DESeq2 normalisation and PSMA trend classification...")

meta_paired$patient_id_f <- factor(meta_paired$patient_id)
meta_paired$type         <- factor(meta_paired$type, levels = c("Initial", "Recurrent"))

dds <- DESeqDataSetFromMatrix(
  countData = counts_paired,
  colData   = meta_paired,
  design    = ~ patient_id_f + type
)
# 'poscounts' handles small sample sizes and potential zero geometric means
dds     <- estimateSizeFactors(dds, type = "poscounts")
vst_mat <- assay(vst(dds, blind = FALSE))

# Store FOLH1 VST values in metadata
meta_paired$vst_val <- as.numeric(vst_mat[folh1_id, meta_paired$sample_barcode])

# Classify each patient by their FOLH1 expression trend at recurrence
dynamics <- meta_paired %>%
  group_by(patient_id) %>%
  summarise(
    vst_ini = vst_val[type == "Initial"],
    vst_rec = vst_val[type == "Recurrent"],
    delta   = vst_rec - vst_ini,
    trend   = ifelse(delta > 0, "Increased", "Decreased"),
    .groups = "drop"
  )

meta_paired <- meta_paired %>%
  left_join(dynamics %>% select(patient_id, trend), by = "patient_id")

n_increased <- sum(dynamics$trend == "Increased")
n_decreased <- sum(dynamics$trend == "Decreased")
message(sprintf("  PSMA Increased: %d patients | PSMA Decreased: %d patients",
                n_increased, n_decreased))

# -----------------------------------------------------------------------------
# 6. Figure 01 — PSMA overview (proportions + paired evolution)
# -----------------------------------------------------------------------------
message("Figure 01: PSMA Overview + Paired T-test...")

p_prop <- ggplot(dynamics, aes(x = trend, fill = trend)) +
  geom_bar(width = 0.6) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  scale_fill_manual(values = c("Increased" = "#D6604D", "Decreased" = "#4393C3")) +
  theme_classic(base_size = 13) +
  labs(title = sprintf("Cohort distribution (n=%d)", nrow(dynamics)),
       x = "PSMA trend", y = "Patient count") +
  theme(legend.position = "none")

p_paired <- ggplot(meta_paired, aes(x = type, y = vst_val)) +
  geom_line(aes(group = patient_id), color = "grey70", alpha = 0.7) +
  geom_boxplot(aes(fill = type), outlier.shape = NA, width = 0.4, alpha = 0.3) +
  geom_point(aes(color = type), size = 3) +
  stat_compare_means(paired = TRUE, method = "t.test", label.x = 1.3) +
  scale_color_manual(values = c("Initial" = "#4393C3", "Recurrent" = "#D6604D")) +
  scale_fill_manual(values  = c("Initial" = "#4393C3", "Recurrent" = "#D6604D")) +
  theme_classic(base_size = 13) +
  labs(title = "Intra-patient FOLH1 evolution",
       y = "FOLH1 VST expression", x = "") +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "01_PSMA_Evolution_Overview.png"),
       p_prop + p_paired, width = 11, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 7. Figure 02 — Volcano plot (DGE at recurrence: Increased vs Decreased)
# -----------------------------------------------------------------------------
message("Figure 02: Differential gene expression at recurrence...")

meta_rec <- meta_paired %>% filter(type == "Recurrent")

if (length(unique(meta_rec$trend)) < 2) {
  message("  Skipping DGE volcano: all recurrent samples have the same trend.")
} else {
  dds_rec <- DESeqDataSetFromMatrix(
    countData = counts_paired[, meta_rec$sample_barcode],
    colData   = meta_rec,
    design    = ~ trend
  )
  dds_rec <- DESeq(dds_rec)
  res_rec <- as.data.frame(
    results(dds_rec, contrast = c("trend", "Increased", "Decreased"))
  ) %>%
    rownames_to_column("gene") %>%
    mutate(sig = ifelse(!is.na(padj) & padj < 0.05, "SIG", "NS"))

  write.csv(res_rec,
            file.path(output_dir, "DGE_Full_Results.csv"),
            row.names = FALSE)

  top_sig <- head(res_rec[res_rec$sig == "SIG", ], 15)
  p_volcano <- ggplot(res_rec, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
    geom_point(alpha = 0.4, size = 1.2) +
    scale_color_manual(values = c("NS" = "grey70", "SIG" = "red3")) +
    geom_text_repel(data = top_sig, aes(label = gene), size = 3, max.overlaps = 30) +
    theme_minimal(base_size = 13) +
    labs(title = "DGE: PSMA Increased vs Decreased (recurrent stage)",
         x = "Log2 fold change", y = "-log10 p-value", color = "")

  ggsave(file.path(output_dir, "02_Volcano_DGE.png"),
         p_volcano, width = 8, height = 6, dpi = 300)
}

# -----------------------------------------------------------------------------
# 8. Figure 03 — Heatmap of top 50 differentially expressed genes
# -----------------------------------------------------------------------------
message("Figure 03: Heatmap of top DEGs...")

if (exists("res_rec") && nrow(res_rec) > 0) {
  top_genes <- head(res_rec$gene[order(res_rec$pvalue)], 50)
  top_genes <- intersect(top_genes, rownames(vst_mat))   # safety check

  if (length(top_genes) >= 2) {
    h_mat <- vst_mat[top_genes, meta_rec$sample_barcode]
    colnames(h_mat) <- meta_rec$patient_id

    png(file.path(output_dir, "03_Heatmap_Top50_DEGs.png"),
        width = 1200, height = 1500, res = 150)
    pheatmap(h_mat,
             annotation_col = data.frame(Trend = meta_rec$trend,
                                         row.names = meta_rec$patient_id),
             main  = "Top 50 DEGs (recurrent samples)",
             scale = "row",
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
    dev.off()
  } else {
    message("  Heatmap skipped: fewer than 2 genes passed the filter.")
  }
}

# -----------------------------------------------------------------------------
# 9. Figure 04 — Biological scores (Δ hypoxia and Δ angiogenesis)
# -----------------------------------------------------------------------------
message("Figure 04: Hallmark pathway scoring...")

hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H")

calc_hallmark_score <- function(vst_data, set_name) {
  ids     <- hallmark_sets %>% filter(gs_name == set_name) %>% pull(ensembl_gene)
  present <- intersect(ids, rownames(vst_data))
  if (length(present) == 0) return(rep(NA_real_, ncol(vst_data)))
  colMeans(t(scale(t(vst_data[present, , drop = FALSE]))), na.rm = TRUE)
}

scores_hyp  <- calc_hallmark_score(vst_mat, "HALLMARK_HYPOXIA")
scores_angio <- calc_hallmark_score(vst_mat, "HALLMARK_ANGIOGENESIS")

meta_paired$hypoxia <- scores_hyp[meta_paired$sample_barcode]
meta_paired$angio   <- scores_angio[meta_paired$sample_barcode]

comp_df <- meta_paired %>%
  group_by(patient_id, trend) %>%
  summarise(
    d_hyp = hypoxia[type == "Recurrent"] - hypoxia[type == "Initial"],
    d_ang = angio[type == "Recurrent"]   - angio[type == "Initial"],
    .groups = "drop"
  )

make_delta_boxplot <- function(data, y_col, title_text) {
  ggplot(data, aes(x = trend, y = .data[[y_col]], fill = trend)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.1, size = 2) +
    stat_compare_means(method = "t.test", label = "p.format", label.x = 1.35) +
    scale_fill_manual(values = c("Increased" = "#D6604D", "Decreased" = "#4393C3")) +
    theme_classic(base_size = 13) +
    labs(title = title_text, x = "PSMA trend", y = "Delta score (Rec - Ini)") +
    theme(legend.position = "none")
}

p_hyp <- make_delta_boxplot(comp_df, "d_hyp", "\u0394 Hypoxia (Rec \u2212 Ini)")
p_ang <- make_delta_boxplot(comp_df, "d_ang", "\u0394 Angiogenesis (Rec \u2212 Ini)")

ggsave(file.path(output_dir, "04_Biological_Scores_TTest.png"),
       p_hyp + p_ang, width = 11, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 10. Figure 05 — Kaplan-Meier overall survival
# -----------------------------------------------------------------------------
message("Figure 05: Kaplan-Meier survival analysis...")

if (!file.exists(clinical_file)) {
  message("  Survival analysis skipped: clinical file not found at:\n  ", clinical_file)
} else {
  clinical <- read.csv(clinical_file) %>%
    mutate(
      patient_id   = bcr_patient_barcode,
      vital_status = ifelse(vital_status == "Dead", 1L, 0L),
      days         = pmax(0,
                          coalesce(suppressWarnings(as.numeric(death_days_to)),
                                   suppressWarnings(as.numeric(last_contact_days_to))),
                          na.rm = TRUE)
    )

  surv_df <- dynamics %>%
    left_join(clinical, by = "patient_id") %>%
    filter(!is.na(days), !is.na(vital_status))

  if (nrow(surv_df) >= 4 && length(unique(surv_df$trend)) == 2) {
    fit <- survfit(Surv(days, vital_status) ~ trend, data = surv_df)

    p_surv <- ggsurvplot(
      fit, data = surv_df,
      pval         = TRUE,
      palette      = c("#4393C3", "#D6604D"),
      title        = "Survival by PSMA expression dynamics",
      legend.labs  = c("PSMA Decreased", "PSMA Increased"),
      xlab         = "Days post-diagnosis",
      ylab         = "Survival probability"
    )
    ggsave(file.path(output_dir, "05_KM_Survival_Analysis.png"),
           p_surv$plot, width = 7, height = 6, dpi = 300)
  } else {
    message("  Survival analysis skipped: insufficient data after merging with clinical file.")
  }
}

# -----------------------------------------------------------------------------
# 11. GO enrichment — genes upregulated in the PSMA Increased group
# -----------------------------------------------------------------------------
message("Step 6: Gene Ontology enrichment (PSMA Increased vs Decreased)...")

if (!exists("res_rec")) {
  message("  GO enrichment skipped: DGE results not available.")
} else {
  up_genes <- res_rec %>%
    filter(log2FoldChange > 0.5, pvalue < 0.05) %>%
    pull(gene)

  if (length(up_genes) < 5) {
    message("  GO enrichment skipped: fewer than 5 upregulated genes detected.")
  } else {
    ego_up <- tryCatch(
      enrichGO(
        gene          = up_genes,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENSEMBL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        readable      = TRUE
      ),
      error = function(e) { message("  GO enrichment error: ", e$message); NULL }
    )

    if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
      write.csv(as.data.frame(ego_up),
                file.path(output_dir, "06_GO_Enrichment_Increased_Group.csv"),
                row.names = FALSE)

      p_dot <- dotplot(ego_up, showCategory = 15) +
        labs(title = "GO terms upregulated in PSMA Increased group")
      p_bar <- barplot(ego_up, showCategory = 15) +
        labs(title = "Top biological processes — PSMA Increased group")

      ggsave(file.path(output_dir, "06_GO_Dotplot.png"),
             p_dot, width = 9, height = 7, dpi = 300)
      ggsave(file.path(output_dir, "07_GO_Barplot.png"),
             p_bar, width = 9, height = 7, dpi = 300)

      message("  GO plots saved.")
    } else {
      message("  No significant GO terms found at the specified thresholds.")
    }
  }
}

# -----------------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------------
message("\n=== DONE. All results saved to: ", output_dir, " ===")
message("  Paired patients analysed: ", length(paired_pts))
message("  PSMA Increased: ", n_increased, " | PSMA Decreased: ", n_decreased)
