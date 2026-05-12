#!/usr/bin/env Rscript

# =============================================================================
# TCGA GBM/LGG - Automated Data Preparation Pipeline (Portable Version)
# =============================================================================

# =============================================================================
# SECTION 0 - CONFIGURATION & PORTABILITY
# =============================================================================

BASE_DIR     <- getwd()
LAVENYA_XLSX <- file.path(BASE_DIR, "11060_2026_5470_MOESM3_ESM.xlsx")

if (!file.exists(LAVENYA_XLSX)) {
  stop(paste("\n[CRITICAL ERROR]: Master Excel file not found at:", LAVENYA_XLSX))
}

OUTPUT_DIR   <- file.path(BASE_DIR, "TCGA_GBM_IDHwt_with_RNAseq")
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

CLINICAL_CSV <- file.path(OUTPUT_DIR, "TCGA_GBM_clinical.csv")

GENES_OF_INTEREST <- c(
  "ENSG00000086205", "ENSG00000112715", "ENSG00000150630", "ENSG00000107159"
  )

cat("✔ Configuration initialized.\n")

# =============================================================================
# SECTION 1 - INSTALL & LOAD PACKAGES
# =============================================================================

cat("=== [1/6] Loading Packages ===\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("TCGAbiolinks", "SummarizedExperiment")
cran_pkgs <- c("readxl", "readr", "dplyr", "tibble")

for (pkg in bioc_pkgs) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE)
for (pkg in cran_pkgs) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr)
  library(tibble); library(SummarizedExperiment); library(TCGAbiolinks)
})

# =============================================================================
# HELPERS
# =============================================================================

extract_patient_id <- function(barcode) substr(trimws(as.character(barcode)), 1, 12)
extract_sample_type_code <- function(barcode) ifelse(nchar(barcode) >= 15, substr(barcode, 14, 15), NA_character_)

sample_type_label <- function(code) {
  dplyr::case_when(code == "01" ~ "Primary Tumor", code == "02" ~ "Recurrent Tumor", TRUE ~ "Other")
}

prepare_expression_df <- function(mat) {
  df <- as.data.frame(t(mat))
  df$sample_barcode <- rownames(df)
  df$patient_id <- extract_patient_id(df$sample_barcode)
  df$sample_type_code <- extract_sample_type_code(df$sample_barcode)
  df$sample_type_label <- sample_type_label(df$sample_type_code)
  df
}

select_gene_columns <- function(df, genes_of_interest) {
  gene_pattern <- paste0("^(", paste(genes_of_interest, collapse = "|"), ")")
  cols_selected <- grep(gene_pattern, colnames(df), value = TRUE)
  if (length(cols_selected) == 0) stop("No genes found.")
  cols_selected
}

write_clinical_subset <- function(clinical_df, keep_patient_ids, out_file, label) {
  # In Biotab format, the column is usually 'bcr_patient_barcode'
  id_col <- if("bcr_patient_barcode" %in% colnames(clinical_df)) "bcr_patient_barcode" else "case_id"
  
  clinical_sub <- clinical_df %>%
    mutate(temp_id_12 = substr(!!sym(id_col), 1, 12)) %>%
    filter(temp_id_12 %in% keep_patient_ids) %>%
    select(-temp_id_12)
  write_csv(clinical_sub, out_file)
  cat(sprintf("  ✔ %s generated.\n", label))
}

write_expression_subset <- function(expr_df, keep_codes, valid_patient_ids, out_file, label, selected_cols = NULL) {
  sub_df <- expr_df %>% filter(sample_type_code %in% keep_codes, patient_id %in% valid_patient_ids)
  if (!is.null(selected_cols)) {
    sub_df <- sub_df[, c("sample_barcode", selected_cols), drop = FALSE]
  } else {
    sub_df <- sub_df[, setdiff(colnames(sub_df), c("patient_id", "sample_type_code", "sample_type_label")), drop = FALSE]
  }
  rownames(sub_df) <- sub_df$sample_barcode
  sub_df$sample_barcode <- NULL
  write.csv(sub_df, file = out_file, row.names = TRUE)
}

# =============================================================================
# SECTION 2 - CLINICAL DATA ACQUISITION (Biotab Method)
# =============================================================================

cat("=== [2/6] Clinical Data Acquisition (GBM & LGG) ===\n")
if (!file.exists(CLINICAL_CSV)) {
  cat("  Downloading Biotab clinical data...\n")
  query_clin <- GDCquery(
    project = c("TCGA-GBM", "TCGA-LGG"), 
    data.category = "Clinical", 
    data.type = "Clinical Supplement", 
    data.format = "BCR Biotab"
  )
  GDCdownload(query_clin)
  clinical_data_list <- GDCprepare(query_clin)
  
  # Extracting only the main patient table from the list
  target_idx <- grep("clinical_patient", names(clinical_data_list))
  all_idx    <- grep("clinical_patient_all", names(clinical_data_list))
  final_idx  <- if(length(all_idx) > 0) all_idx[1] else target_idx[1]
  
  clinical_raw <- clinical_data_list[[final_idx]]
  
  # Security: ensure it's a data frame before writing
  if(is.data.frame(clinical_raw)) {
    write_csv(clinical_raw, CLINICAL_CSV)
    cat("  ✔ Clinical file saved successfully.\n")
  } else {
    stop("Extracted clinical object is not a data frame.")
  }
} else {
  clinical_raw <- read_csv(CLINICAL_CSV, show_col_types = FALSE)
  cat("  ✔ Local clinical file loaded.\n")
}

# =============================================================================
# SECTION 3 - FILTER CLINICAL DATA
# =============================================================================

cat("=== [3/6] Filtering with Master File ===\n")
lavenya_ids_12 <- extract_patient_id(read_excel(LAVENYA_XLSX)$track_name)

# Detecting ID column (Biotab usually uses bcr_patient_barcode)
id_col <- if("bcr_patient_barcode" %in% colnames(clinical_raw)) "bcr_patient_barcode" else "case_id"

clinical_filtered <- clinical_raw %>% 
  filter(substr(!!sym(id_col), 1, 12) %in% lavenya_ids_12)

valid_ids <- unique(substr(clinical_filtered[[id_col]], 1, 12))
out_clinical <- file.path(OUTPUT_DIR, "TCGA_GBM_clinical_inner_join.csv")
write_csv(clinical_filtered, out_clinical)

# =============================================================================
# SECTION 4 - RNA-SEQ DOWNLOAD (Folders & TSV)
# =============================================================================

cat("=== [4/6] Downloading Transcriptome Profiling (STAR - Counts) ===\n")
query_rna <- GDCquery(
  project = c("TCGA-GBM", "TCGA-LGG"), 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

# Downloads folders for each patient containing the .tsv files
GDCdownload(query_rna)
glioma_se <- GDCprepare(query_rna)

# =============================================================================
# SECTION 5 - PROCESSING MATRICES
# =============================================================================

cat("=== [5/6] Processing Matrices ===\n")
# FPKM
fpkm_df <- prepare_expression_df(assay(glioma_se, "fpkm_unstrand"))
fpkm_cols <- select_gene_columns(fpkm_df, GENES_OF_INTEREST)

out_fpkm_primary   <- file.path(OUTPUT_DIR, "fpkm_final_selected_genes.csv")
out_fpkm_all       <- file.path(OUTPUT_DIR, "fpkm_initial_recurrent_selected_genes.csv")
out_fpkm_recurrent <- file.path(OUTPUT_DIR, "fpkm_recurrent_selected_genes.csv")

write_expression_subset(fpkm_df, "01", valid_ids, out_fpkm_primary, "FPKM Pri", fpkm_cols)
write_expression_subset(fpkm_df, c("01", "02"), valid_ids, out_fpkm_all, "FPKM All", fpkm_cols)
write_expression_subset(fpkm_df, "02", valid_ids, out_fpkm_recurrent, "FPKM Rec", fpkm_cols)

# Counts
counts_df <- prepare_expression_df(assay(glioma_se, "unstranded"))
out_counts_primary   <- file.path(OUTPUT_DIR, "rna_seq_raw_counts_filtered.csv")
out_counts_all       <- file.path(OUTPUT_DIR, "rna_seq_raw_counts_initial_recurrent_filtered.csv")
out_counts_recurrent <- file.path(OUTPUT_DIR, "rna_seq_raw_counts_recurrent_filtered.csv")

write_expression_subset(counts_df, "01", valid_ids, out_counts_primary, "Counts Pri")
write_expression_subset(counts_df, c("01", "02"), valid_ids, out_counts_all, "Counts All")
write_expression_subset(counts_df, "02", valid_ids, out_counts_recurrent, "Counts Rec")

# =============================================================================
# SECTION 6 - FINAL CLINICAL SUBSETS & SUMMARY
# =============================================================================

cat("=== [6/6] Finalizing Clinical Subsets ===\n")
out_clinical_all       <- file.path(OUTPUT_DIR, "TCGA_GBM_clinical_initial_recurrent.csv")
out_clinical_recurrent <- file.path(OUTPUT_DIR, "TCGA_GBM_clinical_recurrent_only.csv")

write_clinical_subset(clinical_filtered, valid_ids, out_clinical, "Clinical Join")
write_clinical_subset(clinical_filtered, unique(counts_df$patient_id[counts_df$sample_type_code %in% c("01", "02")]), out_clinical_all, "Clinical All")
write_clinical_subset(clinical_filtered, unique(counts_df$patient_id[counts_df$sample_type_code == "02"]), out_clinical_recurrent, "Clinical Rec")

cat("\n=================================================================\n")
cat("        PIPELINE COMPLETED: TSV DOWNLOADED & BIOTAB FILTERED\n")
cat("=================================================================\n")