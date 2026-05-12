#!/usr/bin/env Rscript

# =============================================================================
# PSMA / FOLH1 analysis in TCGA GBM
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 0. Configurable paths
# -----------------------------------------------------------------------------
# All paths are resolved relative to the script's own directory so the project
# runs on any machine without editing hardcoded paths.
#
# Expected directory layout (created automatically for outputs):
#
#   <project_root>/                              ← where this script lives
#     TCGA_GBM_PSMA_analysis_final.R
#     data/
#       TCGA_GBM_IDHwt_with_RNAseq/
#         TCGA_GBM_clinical_initial_recurrent.csv
#         rna_seq_raw_counts_initial_recurrent_filtered.csv
#         fpkm_initial_recurrent_selected_genes.csv
#     outputs/
#       PSMA_GBM_initial_recurrent_outputs/      ← generated automatically
#         Initial/
#         Recurrent/
#         Paired_Analysis/
#         Initial_vs_Recurrent/
#         PSM_1to2/
#
# ── Detect the script's own directory ────────────────────────────────────────
get_script_dir <- function() {
  # RStudio interactive / sourced file
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  # Rscript --file= (command-line execution)
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0)
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  # Fallback: current working directory
  getwd()
}

script_dir <- get_script_dir()
message("Script directory: ", script_dir)

# ── Input data ────────────────────────────────────────────────────────────────
data_dir      <- file.path(script_dir, "data", "TCGA_GBM_IDHwt_with_RNAseq")
clinical_file <- file.path(data_dir, "TCGA_GBM_clinical_initial_recurrent.csv")
counts_file   <- file.path(data_dir, "rna_seq_raw_counts_initial_recurrent_filtered.csv")
fpkm_file     <- file.path(data_dir, "fpkm_initial_recurrent_selected_genes.csv")

# ── Output root (sub-folders created automatically) ───────────────────────────
output_root   <- file.path(script_dir, "outputs", "PSMA_GBM_initial_recurrent_outputs")

folh1_id  <- "ENSG00000086205"
psm_ratio <- 2L   # 2 initial matched per recurrent (1:2 PSM)

# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("DESeq2", "apeglm", "org.Hs.eg.db", "AnnotationDbi",
               "clusterProfiler", "enrichplot")
cran_pkgs <- c("dplyr", "tibble", "tidyr", "ggplot2", "survival",
               "gprofiler2", "survminer", "msigdbr", "scales", "ggrepel",
               "MatchIt", "patchwork")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

suppressPackageStartupMessages({
  library(tibble); library(tidyr); library(ggplot2); library(gprofiler2)
  library(survival); library(survminer); library(msigdbr); library(scales)
  library(ggrepel); library(org.Hs.eg.db); library(AnnotationDbi)
  library(clusterProfiler); library(enrichplot); library(DESeq2)
  library(dplyr); library(MatchIt); library(patchwork)
})

select <- dplyr::select; filter <- dplyr::filter
rename <- dplyr::rename; mutate <- dplyr::mutate

# -----------------------------------------------------------------------------
# 2. General helpers
# -----------------------------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

clean_ensembl <- function(x) gsub("\\.\\d+$", "", x)

sample_type_label <- function(code) {
  dplyr::case_when(code == "01" ~ "Initial", code == "02" ~ "Recurrent",
                   TRUE ~ paste0("Other_", code))
}

prepare_sample_metadata <- function(barcodes) {
  data.frame(
    sample_barcode   = barcodes,
    patient_id       = substr(barcodes, 1, 12),
    sample_type_code = substr(barcodes, 14, 15),
    sample_type      = sample_type_label(substr(barcodes, 14, 15)),
    sample_key       = paste(substr(barcodes, 1, 12),
                             substr(barcodes, 14, 15), sep = "_"),
    stringsAsFactors = FALSE
  )
}

deduplicate_by_patient_and_type <- function(df, id_vector) {
  meta     <- prepare_sample_metadata(id_vector)
  keep_idx <- !duplicated(meta$sample_key)
  list(data      = df[keep_idx, , drop = FALSE],
       meta      = meta[keep_idx, , drop = FALSE],
       removed_n = sum(!keep_idx))
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

# Converts TCGA "[Not Available]" / "[Not Applicable]" values to NA
clean_tcga_val <- function(x) {
  ifelse(grepl("\\[Not", as.character(x)), NA_character_, as.character(x))
}

safe_wilcox <- function(formula, data) {
  tryCatch(wilcox.test(formula, data = data), error = function(e) NULL)
}

safe_wilcox_groups <- function(values, groups, group_a, group_b) {
  idx    <- !is.na(values) & !is.na(groups) & groups %in% c(group_a, group_b)
  vals_a <- values[idx][as.character(groups[idx]) == group_a]
  vals_b <- values[idx][as.character(groups[idx]) == group_b]
  if (length(vals_a) == 0 || length(vals_b) == 0) return(NULL)
  tryCatch(wilcox.test(vals_a, vals_b), error = function(e) NULL)
}

# T-test between two groups - returns a formatted p-value string for plot annotation
safe_ttest_groups <- function(values, groups, group_a, group_b) {
  idx    <- !is.na(values) & !is.na(groups) & groups %in% c(group_a, group_b)
  vals_a <- values[idx][as.character(groups[idx]) == group_a]
  vals_b <- values[idx][as.character(groups[idx]) == group_b]
  if (length(vals_a) < 2 || length(vals_b) < 2) return(NA_real_)
  tryCatch(t.test(vals_a, vals_b)$p.value, error = function(e) NA_real_)
}

format_pval <- function(p) {
  if (is.na(p))      return("p = NA")
  if (p < 0.001)     return(sprintf("p = %.2e", p))
  if (p < 0.05)      return(sprintf("p = %.3f", p))
  return(sprintf("p = %.3f", p))
}

# Adds a t-test annotation to a ggplot boxplot comparing two groups
add_ttest_annotation <- function(p, df, y_col, group_col, group_a, group_b,
                                 y_offset_frac = 0.08) {
  pval <- safe_ttest_groups(df[[y_col]], df[[group_col]], group_a, group_b)
  y_max <- max(df[[y_col]], na.rm = TRUE)
  y_min <- min(df[[y_col]], na.rm = TRUE)
  y_ann <- y_max + (y_max - y_min) * y_offset_frac
  p + annotate("text", x = 1.5, y = y_ann,
               label = format_pval(pval),
               size = 4, fontface = "italic", color = "black")
}

# Universal volcano helper - replaces all repetitive volcano plotting code
# Gère les padj=0 (Inf), adapte l'axe Y automatiquement
row_zscore <- function(mat) {
  z <- t(scale(t(as.matrix(mat))))
  z[is.nan(z)] <- 0
  z[is.infinite(z)] <- NA_real_
  z
}

cluster_order_from_z <- function(z_mat, margin = c("rows", "cols")) {
  margin <- match.arg(margin)
  x <- if (margin == "rows") as.matrix(z_mat) else t(as.matrix(z_mat))
  if (nrow(x) < 2) return(rownames(x))
  
  x[!is.finite(x)] <- NA_real_
  row_means <- rowMeans(x, na.rm = TRUE)
  row_means[!is.finite(row_means)] <- 0
  for (i in seq_len(nrow(x))) {
    x[i, is.na(x[i, ])] <- row_means[i]
  }
  
  hc <- tryCatch(hclust(dist(x)), error = function(e) NULL)
  if (is.null(hc)) rownames(x) else rownames(x)[hc$order]
}

z_matrix_to_long <- function(z_mat, sample_order = NULL, gene_order = NULL) {
  z_mat <- as.matrix(z_mat)
  if (is.null(sample_order))
    sample_order <- cluster_order_from_z(z_mat, "cols")
  if (is.null(gene_order))
    gene_order <- cluster_order_from_z(z_mat[, sample_order, drop = FALSE], "rows")
  
  z_sub <- z_mat[gene_order, sample_order, drop = FALSE]
  data.frame(gene_label = rownames(z_sub), z_sub, check.names = FALSE) %>%
    tidyr::pivot_longer(-gene_label, names_to = "sample_barcode", values_to = "z") %>%
    mutate(sample_barcode = factor(sample_barcode, levels = sample_order),
           gene_label     = factor(gene_label, levels = rev(gene_order)))
}

z_heatmap_fill <- function(name = "z-score") {
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = name, na.value = "grey90")
}

make_volcano <- function(df,
                         lfc_col       = "log2FoldChange",
                         padj_col      = "padj",
                         label_col     = "SYMBOL",
                         gene_id_col   = "gene_id",
                         lfc_threshold = 1,
                         padj_threshold = 0.05,
                         up_label      = "Up",
                         down_label    = "Down",
                         n_labels      = 20,
                         title         = "",
                         subtitle      = "",
                         x_label       = "log2FC (apeglm)",
                         out_file      = NULL,
                         width = 7, height = 5) {
  
  vdf <- df %>%
    filter(!is.na(.data[[padj_col]]), !is.na(.data[[lfc_col]])) %>%
    mutate(
      label_text = dplyr::if_else(
        is.na(.data[[label_col]]) | .data[[label_col]] == "",
        .data[[gene_id_col]], .data[[label_col]]),
      # pmax prevents Inf when padj == 0 due to floating-point underflow
      nlp = -log10(pmax(.data[[padj_col]], 1e-300)),
      sig = dplyr::case_when(
        .data[[padj_col]] < padj_threshold & .data[[lfc_col]] >  lfc_threshold ~ up_label,
        .data[[padj_col]] < padj_threshold & .data[[lfc_col]] < -lfc_threshold ~ down_label,
        TRUE ~ "NS"
      )
    )
  
  # Axis Y : clip at 95th pct of significant genes × 1.35 (at least 2× sig line)
  sig_nlp <- vdf$nlp[vdf$sig != "NS" & is.finite(vdf$nlp)]
  y_cap   <- if (length(sig_nlp) > 0)
    max(quantile(sig_nlp, 0.95, na.rm = TRUE) * 1.35,
        -log10(padj_threshold) * 2.5)
  else
    max(vdf$nlp[is.finite(vdf$nlp)], na.rm = TRUE) * 1.1
  
  n_beyond <- sum(vdf$nlp > y_cap & is.finite(vdf$nlp), na.rm = TRUE)
  if (n_beyond > 0) {
    note <- sprintf("%d gene(s) not shown (above y-axis limit; shown in _DESeq2_all.csv)",
                    n_beyond)
    subtitle <- if (nchar(subtitle) > 0) paste0(subtitle, "\n", note) else note
  }
  
  # Labels: only genes within the visible area, top N by padj
  top_df <- vdf %>%
    filter(sig != "NS", nlp <= y_cap * 1.02) %>%
    arrange(.data[[padj_col]]) %>%
    slice_head(n = n_labels)
  
  col_vals <- setNames(c("#D6604D", "#4393C3", "grey70"),
                       c(up_label, down_label, "NS"))
  # Drop colour levels not present (prevents legend entries for empty groups)
  present   <- intersect(c(up_label, down_label, "NS"), unique(vdf$sig))
  col_vals  <- col_vals[present]
  
  p <- ggplot(vdf, aes(x = .data[[lfc_col]], y = nlp, color = sig)) +
    geom_point(alpha = 0.45, size = 1.2) +
    geom_label_repel(data = top_df, aes(label = label_text),
                     size = 2.5, max.overlaps = 30, show.legend = FALSE,
                     box.padding = 0.3, segment.size = 0.3) +
    scale_color_manual(values = col_vals, drop = FALSE) +
    coord_cartesian(ylim = c(0, y_cap * 1.05)) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = "dashed", linewidth = 0.4, color = "grey40") +
    geom_hline(yintercept = -log10(padj_threshold),
               linetype = "dashed", linewidth = 0.4, color = "grey40") +
    theme_classic(base_size = 13) +
    theme(plot.subtitle = element_text(size = 8, color = "grey40")) +
    labs(title = title, subtitle = subtitle,
         x = x_label, y = "-log10 adj. p", color = "")
  
  if (!is.null(out_file))
    save_plot(out_file, p, width = width, height = height)
  
  invisible(p)
}

# Actual clinical columns in the TCGA GBM file:
#   bcr_patient_barcode, death_days_to, last_contact_days_to,
#   vital_status, age_at_initial_pathologic_diagnosis,
#   gender (uppercase), pharmaceutical_tx_adjuvant, radiation_treatment_adjuvant
prepare_survival_df <- function(clinical_df) {
  clinical_df %>%
    mutate(
      vital_status    = as.character(vital_status),
      death_event     = as.integer(vital_status == "Dead"),
      death_days_num  = safe_numeric(clean_tcga_val(death_days_to)),
      follow_days_num = safe_numeric(clean_tcga_val(last_contact_days_to)),
      OS_days = dplyr::case_when(
        !is.na(death_days_num)  & death_days_num  > 0 ~ death_days_num,
        !is.na(follow_days_num) & follow_days_num > 0 ~ follow_days_num,
        TRUE ~ NA_real_),
      OS_months = OS_days / 30.44
    ) %>%
    select(bcr_patient_barcode, vital_status, death_event, OS_days, OS_months) %>%
    rename(case_id = bcr_patient_barcode) %>%
    distinct(case_id, .keep_all = TRUE)
}

make_psma_groups <- function(expr_values) {
  if (length(expr_values) < 2 || length(unique(expr_values)) < 2)
    stop("FOLH1 expression has insufficient variation to define PSMA groups.")
  ranks  <- rank(expr_values, ties.method = "first")
  cutoff <- floor(length(expr_values) / 2)
  factor(ifelse(ranks <= cutoff, "PSMA Low", "PSMA High"),
         levels = c("PSMA Low", "PSMA High"))
}

map_ensembl_to_annotations <- function(ensembl_ids) {
  AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(ensembl_ids),
    columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = "ENSEMBL"
  ) %>% distinct(ENSEMBL, .keep_all = TRUE)
}

save_plot <- function(filename, plot_obj, width = 6, height = 4.5, dpi = 300) {
  ggsave(filename = filename, plot = plot_obj,
         width = width, height = height, dpi = dpi)
}

write_gene_matrix <- function(mat, gene_map, out_file) {
  as.data.frame(mat) %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_map, by = c("gene_id" = "ENSEMBL")) %>%
    relocate(gene_id, SYMBOL) %>%
    write.csv(out_file, row.names = FALSE)
}

make_survival_formula <- function(group_var) {
  as.formula(sprintf("Surv(OS_days, death_event) ~ %s", group_var))
}

get_coef_name <- function(dds, pattern) {
  hit <- grep(pattern, resultsNames(dds), value = TRUE)
  if (length(hit) == 0) stop(sprintf("No DESeq2 coef matched: %s", pattern))
  hit[1]
}

# -----------------------------------------------------------------------------
# 3. Helpers for GO plots
# -----------------------------------------------------------------------------

flatten_list_cols <- function(df) {
  list_cols <- vapply(df, is.list, logical(1))
  df[list_cols] <- lapply(df[list_cols], function(col) {
    vapply(col, function(x) {
      if (is.null(x) || length(x) == 0) NA_character_
      else paste(as.character(x), collapse = ",")
    }, character(1))
  })
  df
}

safe_enrich_barplot <- function(enrich_res, n_show = 15, title = "") {
  res_df <- as.data.frame(enrich_res)
  if (nrow(res_df) == 0) return(invisible(NULL))
  plot_df <- head(res_df[order(res_df$p.adjust), ], n_show)
  plot_df$Description <- factor(plot_df$Description,
                                levels = rev(plot_df$Description))
  if (!"Count" %in% colnames(plot_df))
    plot_df$Count <- as.integer(sub("/.*", "", plot_df$GeneRatio))
  ggplot(plot_df, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "#D6604D", high = "#4393C3", name = "adj. p",
                        guide = guide_colorbar(reverse = TRUE)) +
    theme_classic(base_size = 11) +
    theme(axis.text.y   = element_text(size = 9),
          plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8, color = "grey40")) +
    labs(title = title, x = "Gene count", y = NULL)
}

safe_enrich_dotplot <- function(enrich_res, n_show = 15, title = "") {
  res_df <- as.data.frame(enrich_res)
  if (nrow(res_df) == 0) return(invisible(NULL))
  plot_df <- head(res_df[order(res_df$p.adjust), ], n_show)
  plot_df$Description <- factor(plot_df$Description,
                                levels = rev(plot_df$Description))
  if (!"Count" %in% colnames(plot_df))
    plot_df$Count <- as.integer(sub("/.*", "", plot_df$GeneRatio))
  plot_df$GeneRatio_num <- plot_df$Count /
    as.integer(sub(".*/", "", plot_df$GeneRatio))
  ggplot(plot_df, aes(x = GeneRatio_num, y = Description,
                      color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "#D6604D", high = "#4393C3", name = "adj. p",
                         guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 9),
          plot.title  = element_text(size = 10, face = "bold")) +
    labs(title = title, x = "Gene ratio", y = NULL)
}

# -----------------------------------------------------------------------------
# 4. run_go_enrichment
# -----------------------------------------------------------------------------

run_go_enrichment <- function(sig_df, gene_map, cohort_dir, prefix,
                              positive_label, negative_label,
                              logfc_col = "log2FoldChange", padj_col = "padj",
                              gene_id_col = "gene_id", logfc_threshold = 0.5,
                              min_genes = 10) {
  
  if (is.null(sig_df) || nrow(sig_df) == 0) {
    writeLines("GO skipped: no significant genes.",
               file.path(cohort_dir, paste0(prefix, "_GO_summary.txt")))
    return(invisible(NULL))
  }
  
  go_input <- sig_df %>%
    filter(!is.na(.data[[padj_col]]), .data[[padj_col]] < 0.05,
           !is.na(.data[[logfc_col]]))
  
  if (nrow(go_input) == 0) {
    writeLines("GO skipped: no genes passed padj < 0.05.",
               file.path(cohort_dir, paste0(prefix, "_GO_summary.txt")))
    return(invisible(NULL))
  }
  
  all_sig_genes <- go_input %>% pull(.data[[gene_id_col]]) %>% unique()
  genes_up      <- go_input %>%
    filter(.data[[logfc_col]] >  logfc_threshold) %>%
    pull(.data[[gene_id_col]]) %>% unique()
  genes_down    <- go_input %>%
    filter(.data[[logfc_col]] < -logfc_threshold) %>%
    pull(.data[[gene_id_col]]) %>% unique()
  
  writeLines(
    c(paste("Prefix:", prefix),
      paste("All sig DEGs:", length(all_sig_genes)),
      paste(sprintf("Up (logFC > %g):", logfc_threshold), length(genes_up)),
      paste(sprintf("Down (logFC < -%g):", logfc_threshold), length(genes_down))),
    file.path(cohort_dir, paste0(prefix, "_GO_summary.txt")))
  
  run_full_go <- function(gene_list, group_name) {
    if (length(gene_list) < min_genes) {
      writeLines(sprintf("GO skipped for '%s': %d genes (min %d).",
                         group_name, length(gene_list), min_genes),
                 file.path(cohort_dir, paste0(prefix, "_GO_", group_name, "_skipped.txt")))
      return(invisible(NULL))
    }
    ont_labels <- c(ALL = "All ontologies",   BP = "Biological Process",
                    MF  = "Molecular Function", CC = "Cellular Component")
    go_list <- list()
    for (ont in names(ont_labels)) {
      res <- tryCatch(
        enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                 ont = ont, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                 readable = FALSE),
        error = function(e) NULL)
      if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
        go_list[[ont]] <- res
        write.csv(as.data.frame(res),
                  file.path(cohort_dir,
                            paste0(prefix, "_GO_", ont, "_", group_name, "_terms.csv")),
                  row.names = FALSE)
      }
    }
    if (length(go_list) == 0) {
      writeLines(sprintf("No GO enrichment found for '%s'.", group_name),
                 file.path(cohort_dir, paste0(prefix, "_GO_", group_name, "_no_results.txt")))
      return(invisible(NULL))
    }
    pdf_path <- file.path(cohort_dir,
                          paste0(prefix, "_GO_barplots_", group_name, ".pdf"))
    pdf(pdf_path, width = 11, height = 8)
    for (ont in names(ont_labels)) {
      if (!ont %in% names(go_list)) next
      res        <- go_list[[ont]]
      n_show     <- min(15, nrow(as.data.frame(res)))
      title_base <- sprintf("GO %s — %s\n[%s | %d genes | %d terms]",
                            ont_labels[ont], group_name, prefix,
                            length(gene_list), nrow(as.data.frame(res)))
      p_bar <- safe_enrich_barplot(res, n_show, title_base)
      if (!is.null(p_bar)) print(p_bar)
      p_dot <- safe_enrich_dotplot(res, n_show, paste("Dotplot |", title_base))
      if (!is.null(p_dot)) print(p_dot)
    }
    dev.off()
    message(sprintf("  GO PDF saved -> %s", pdf_path))
  }
  
  run_gprofiler_analysis <- function(gene_list, group_name) {
    if (length(gene_list) <= 5) return(invisible(NULL))
    gost_res <- tryCatch(
      gost(query = gene_list, organism = "hsapiens", significant = TRUE,
           exclude_iea = TRUE, user_threshold = 0.05,
           correction_method = "g_SCS",
           sources = c("GO:BP", "KEGG", "REAC")),
      error = function(e) NULL)
    if (is.null(gost_res) || is.null(gost_res$result) ||
        nrow(gost_res$result) == 0) {
      writeLines(sprintf("No g:Profiler enrichment for '%s'.", group_name),
                 file.path(cohort_dir,
                           paste0(prefix, "_gProfiler_", group_name, "_summary.txt")))
      return(invisible(NULL))
    }
    write.csv(flatten_list_cols(gost_res$result),
              file.path(cohort_dir,
                        paste0(prefix, "_gProfiler_Results_", group_name, ".csv")),
              row.names = FALSE)
    tryCatch({
      p      <- gostplot(gost_res, interactive = FALSE, capped = TRUE)
      p_labs <- publish_gostplot(
        p, highlight_terms = head(unique(gost_res$result$term_id), 10))
      ggsave(file.path(cohort_dir,
                       paste0(prefix, "_gProfiler_Manhattan_", group_name, ".pdf")),
             p_labs, width = 12, height = 8)
    }, error = function(e) {
      message(sprintf("  gProfiler plot failed for '%s': %s", group_name, e$message))
    })
  }
  
  run_full_go(all_sig_genes, "all_sig")
  run_gprofiler_analysis(all_sig_genes, "all_sig")
  run_full_go(genes_up,   paste0(positive_label, "_High"))
  run_gprofiler_analysis(genes_up, paste0(positive_label, "_High"))
  run_full_go(genes_down, paste0(negative_label, "_High"))
  run_gprofiler_analysis(genes_down, paste0(negative_label, "_High"))
}

# -----------------------------------------------------------------------------
# 5. run_hypoxia_analysis
# -----------------------------------------------------------------------------

run_hypoxia_analysis <- function(vsd_mat, cohort_meta, hallmark_ids,
                                 hallmark_symbols, cohort_dir, cohort_label,
                                 group_var    = "PSMA_group",
                                 group_levels = c("PSMA Low", "PSMA High"),
                                 palette      = c("PSMA Low"  = "#4393C3",
                                                  "PSMA High" = "#D6604D")) {
  
  hypoxia_ids <- intersect(hallmark_ids, rownames(vsd_mat))
  if (length(hypoxia_ids) == 0) {
    writeLines("Hypoxia analysis skipped: no Hallmark genes in VST matrix.",
               file.path(cohort_dir, paste0(cohort_label, "_hypoxia_summary.txt")))
    return(invisible(NULL))
  }
  
  z_mat         <- row_zscore(vsd_mat[hypoxia_ids, , drop = FALSE])
  hypoxia_score <- colMeans(z_mat, na.rm = TRUE)
  cohort_meta$Hallmark_Hypoxia_score <- hypoxia_score[cohort_meta$sample_barcode]
  
  hypoxia_df <- cohort_meta %>%
    select(sample_barcode, patient_id, sample_type,
           all_of(group_var), FOLH1_vst, Hallmark_Hypoxia_score)
  write.csv(hypoxia_df,
            file.path(cohort_dir, paste0(cohort_label, "_hypoxia_scores.csv")),
            row.names = FALSE)
  
  wt    <- safe_wilcox(as.formula(paste("Hallmark_Hypoxia_score ~", group_var)),
                       hypoxia_df)
  grp_a <- group_levels[1]; grp_b <- group_levels[2]
  write.csv(data.frame(
    cohort          = cohort_label,
    n_hypoxia_genes = length(hypoxia_ids),
    group_a  = grp_a,
    median_a = median(hypoxia_df$Hallmark_Hypoxia_score[hypoxia_df[[group_var]] == grp_a],
                      na.rm = TRUE),
    group_b  = grp_b,
    median_b = median(hypoxia_df$Hallmark_Hypoxia_score[hypoxia_df[[group_var]] == grp_b],
                      na.rm = TRUE),
    wilcox_p = if (is.null(wt)) NA_real_ else wt$p.value),
    file.path(cohort_dir, paste0(cohort_label, "_hypoxia_group_stats.csv")),
    row.names = FALSE)
  
  sym_map      <- map_ensembl_to_annotations(hypoxia_ids)
  gene_list_df <- data.frame(
    ENSEMBL = hypoxia_ids,
    SYMBOL  = sym_map$SYMBOL[match(hypoxia_ids, sym_map$ENSEMBL)])
  write.csv(gene_list_df,
            file.path(cohort_dir, paste0(cohort_label, "_hypoxia_gene_list.csv")),
            row.names = FALSE)
  
  gene_syms_measured <- na.omit(gene_list_df$SYMBOL)
  subtitle_text <- paste0(
    "HALLMARK_HYPOXIA — ", length(hypoxia_ids), " genes measured",
    " (of ", length(hallmark_symbols), " in set) | ",
    "Top: ", paste(head(sort(gene_syms_measured), 10), collapse = ", "),
    if (length(gene_syms_measured) > 10) ", ..." else "")
  
  hypoxia_df[[group_var]] <- factor(hypoxia_df[[group_var]], levels = group_levels)
  
  p_box <- ggplot(hypoxia_df,
                  aes(x = .data[[group_var]], y = Hallmark_Hypoxia_score,
                      fill = .data[[group_var]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.55) +
    scale_fill_manual(values = palette) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.subtitle   = element_text(size = 7, color = "grey40"),
          plot.margin     = margin(t = 20, r = 10, b = 10, l = 10)) +
    labs(title    = paste0(cohort_label, " — Hallmark Hypoxia score"),
         subtitle = subtitle_text,
         x = "", y = "Hypoxia score (mean z-score)")
  p_box <- add_ttest_annotation(p_box, hypoxia_df, "Hallmark_Hypoxia_score",
                                group_var, group_levels[1], group_levels[2])
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_hypoxia_boxplot.png")),
            p_box, width = 6.5, height = 5.5)
  
  p_scatter <- ggplot(hypoxia_df,
                      aes(x = FOLH1_vst, y = Hallmark_Hypoxia_score,
                          color = .data[[group_var]])) +
    geom_point(size = 2.1, alpha = 0.8) +
    scale_color_manual(values = palette) +
    theme_classic(base_size = 13) +
    labs(title = paste0(cohort_label, " FOLH1 vs hypoxia score"),
         x = "FOLH1 VST expression", y = "Hallmark hypoxia score", color = "")
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_FOLH1_vs_hypoxia_scatter.png")),
            p_scatter, width = 5.8, height = 4.5)
  
  gene_var  <- apply(z_mat, 1, var, na.rm = TRUE)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(25, length(gene_var))]
  top_z     <- z_mat[top_genes, cohort_meta$sample_barcode, drop = FALSE]
  top_syms  <- map_ensembl_to_annotations(top_genes)
  rownames(top_z) <- dplyr::if_else(
    top_genes %in% top_syms$ENSEMBL,
    top_syms$SYMBOL[match(top_genes, top_syms$ENSEMBL)],
    top_genes)
  rownames(top_z) <- make.unique(ifelse(is.na(rownames(top_z)) | rownames(top_z) == "",
                                        top_genes, rownames(top_z)))
  
  sample_order <- cohort_meta %>%
    filter(sample_barcode %in% colnames(top_z)) %>%
    arrange(FOLH1_vst) %>%
    pull(sample_barcode)
  if (length(sample_order) == 0)
    sample_order <- cluster_order_from_z(top_z, "cols")
  
  marker_order <- cluster_order_from_z(top_z[, sample_order, drop = FALSE], "rows")
  display_z    <- top_z[marker_order, sample_order, drop = FALSE]
  has_folh1_row <- folh1_id %in% rownames(vsd_mat)
  if (has_folh1_row) {
    folh1_z <- row_zscore(vsd_mat[folh1_id, sample_order, drop = FALSE])
    rownames(folh1_z) <- "FOLH1"
    display_z <- rbind(folh1_z, display_z)
  }
  display_order <- rownames(display_z)
  
  heatmap_long <- z_matrix_to_long(display_z, sample_order = sample_order,
                                   gene_order = display_order) %>%
    left_join(cohort_meta %>% select(sample_barcode, all_of(group_var)),
              by = "sample_barcode") %>%
    rename(gene_symbol = gene_label)
  
  p_heat <- ggplot(heatmap_long,
                   aes(x = sample_barcode, y = gene_symbol, fill = z)) +
    geom_tile() +
    z_heatmap_fill() +
    theme_minimal(base_size = 10) +
    theme(axis.text.x   = element_blank(),
          axis.ticks.x  = element_blank(),
          panel.grid    = element_blank(),
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, color = "grey40")) +
    labs(title    = paste0(cohort_label, " — Top 25 most variable Hallmark Hypoxia genes"),
         subtitle = paste0("Ordered by ", group_var,
                           " then hypoxia score | Full list in _hypoxia_gene_list.csv"),
         x = paste0("Samples (n=", ncol(top_z), ")"), y = "Gene")
  p_heat <- p_heat +
    labs(title    = paste0(cohort_label, " - FOLH1 + top 25 Hallmark Hypoxia genes"),
         subtitle = paste0("Samples ordered by FOLH1 expression; hypoxia markers clustered by z-score | ",
                           "Full list in _hypoxia_gene_list.csv"))
  if (has_folh1_row && nrow(display_z) > 1) {
    p_heat <- p_heat +
      geom_hline(yintercept = nrow(display_z) - 0.5,
                 color = "grey25", linewidth = 0.35)
  }
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_hypoxia_top25_heatmap.png")),
            p_heat, width = 9, height = 6)
}

# -----------------------------------------------------------------------------
# 6. run_fpkm_selected_gene_analysis
# -----------------------------------------------------------------------------

run_fpkm_selected_gene_analysis <- function(fpkm_mat, cohort_meta, gene_map,
                                            cohort_dir, cohort_label) {
  if (is.null(fpkm_mat) || nrow(fpkm_mat) == 0) {
    writeLines("FPKM analysis skipped.",
               file.path(cohort_dir, paste0(cohort_label, "_fpkm_summary.txt")))
    return(invisible(NULL))
  }
  common <- intersect(cohort_meta$sample_barcode, rownames(fpkm_mat))
  if (length(common) < 2) {
    writeLines("FPKM analysis skipped: < 2 overlapping samples.",
               file.path(cohort_dir, paste0(cohort_label, "_fpkm_summary.txt")))
    return(invisible(NULL))
  }
  fpkm_sub <- fpkm_mat[common, , drop = FALSE]
  meta_sub  <- cohort_meta %>%
    filter(sample_barcode %in% common) %>%
    arrange(match(sample_barcode, common))
  
  write.csv(as.data.frame(fpkm_sub) %>%
              rownames_to_column("sample_barcode") %>%
              left_join(meta_sub, by = "sample_barcode"),
            file.path(cohort_dir, paste0(cohort_label, "_fpkm_matrix.csv")),
            row.names = FALSE)
  
  long_df <- as.data.frame(fpkm_sub) %>%
    rownames_to_column("sample_barcode") %>%
    tidyr::pivot_longer(-sample_barcode, names_to = "gene_id", values_to = "FPKM") %>%
    left_join(meta_sub, by = "sample_barcode") %>%
    left_join(gene_map, by = c("gene_id" = "ENSEMBL")) %>%
    mutate(gene_label = ifelse(is.na(SYMBOL) | SYMBOL == "", gene_id, SYMBOL),
           log2_FPKM  = log2(FPKM + 1))
  
  stats_df <- long_df %>%
    group_by(gene_id, gene_label) %>%
    summarise(
      median_low  = median(FPKM[PSMA_group == "PSMA Low"],  na.rm = TRUE),
      median_high = median(FPKM[PSMA_group == "PSMA High"], na.rm = TRUE),
      log2FC  = log2(median_high + 1) - log2(median_low + 1),
      p_value = { wt <- safe_wilcox_groups(FPKM, PSMA_group, "PSMA Low", "PSMA High")
      if (is.null(wt)) NA_real_ else wt$p.value },
      .groups = "drop") %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    arrange(padj, desc(abs(log2FC)))
  write.csv(stats_df,
            file.path(cohort_dir, paste0(cohort_label, "_fpkm_stats.csv")),
            row.names = FALSE)
  
  heat_mat <- t(log2(fpkm_sub + 1))
  heat_gene_ids <- rownames(heat_mat)
  heat_labels <- gene_map$SYMBOL[match(heat_gene_ids, gene_map$ENSEMBL)]
  rownames(heat_mat) <- make.unique(ifelse(is.na(heat_labels) | heat_labels == "",
                                           heat_gene_ids, heat_labels))
  heat_z <- row_zscore(heat_mat)
  ordered_s <- cluster_order_from_z(heat_z, "cols")
  ordered_g <- cluster_order_from_z(heat_z[, ordered_s, drop = FALSE], "rows")
  hmap_df <- z_matrix_to_long(heat_z, sample_order = ordered_s,
                              gene_order = ordered_g)
  
  p_heat <- ggplot(hmap_df, aes(sample_barcode, gene_label, fill = z)) +
    geom_tile() +
    z_heatmap_fill() +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid  = element_blank()) +
    labs(title = paste0(cohort_label, " — selected hypoxic genes (FPKM)"),
         subtitle = "Rows and samples clustered by row z-score",
         x = "Samples", y = "Gene", fill = "z-score")
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_fpkm_heatmap.png")),
            p_heat, width = 8, height = 4.8)
}

# -----------------------------------------------------------------------------
# 7. run_deseq_analysis
# -----------------------------------------------------------------------------

run_deseq_analysis <- function(count_mat, sample_meta, gene_map,
                               cohort_dir, cohort_label) {
  grp_counts <- table(sample_meta$PSMA_group)
  if (length(grp_counts) < 2 || any(grp_counts < 2)) {
    writeLines("DESeq2 skipped: need >= 2 samples per PSMA group.",
               file.path(cohort_dir, paste0(cohort_label, "_DESeq2_summary.txt")))
    return(invisible(NULL))
  }
  dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                                colData   = sample_meta,
                                design    = ~ PSMA_group)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  
  res        <- results(dds, contrast = c("PSMA_group", "PSMA High", "PSMA Low"),
                        alpha = 0.05)
  res_shrunk <- lfcShrink(dds, coef = get_coef_name(dds, "^PSMA_group_"),
                          type = "apeglm")
  
  ann <- function(r) as.data.frame(r) %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_map, by = c("gene_id" = "ENSEMBL")) %>%
    arrange(padj)
  
  res_df        <- ann(res)
  res_shrunk_df <- ann(res_shrunk)
  sig_df        <- res_shrunk_df %>% filter(!is.na(padj), padj < 0.05)
  
  write.csv(res_df,
            file.path(cohort_dir, paste0(cohort_label, "_DESeq2_all_genes.csv")),
            row.names = FALSE)
  write.csv(res_shrunk_df,
            file.path(cohort_dir, paste0(cohort_label, "_DESeq2_LFCshrink.csv")),
            row.names = FALSE)
  write.csv(sig_df,
            file.path(cohort_dir, paste0(cohort_label, "_DESeq2_sig_FDR0.05.csv")),
            row.names = FALSE)
  write.csv(sig_df %>% filter(log2FoldChange > 0),
            file.path(cohort_dir, paste0(cohort_label, "_DESeq2_up_PSMA_High.csv")),
            row.names = FALSE)
  write.csv(sig_df %>% filter(log2FoldChange < 0),
            file.path(cohort_dir, paste0(cohort_label, "_DESeq2_down_PSMA_High.csv")),
            row.names = FALSE)
  
  make_volcano(
    df         = res_shrunk_df,
    up_label   = "Up in PSMA High",
    down_label = "Down in PSMA High",
    title      = paste0(cohort_label, " — PSMA High vs Low"),
    out_file   = file.path(cohort_dir, paste0(cohort_label, "_DESeq2_volcano.png"))
  )
  
  writeLines(c(paste("Cohort:", cohort_label),
               paste("Samples:", ncol(count_mat)),
               paste("Genes tested:", nrow(res_df)),
               paste("Sig FDR<0.05:", nrow(sig_df)),
               paste("Up PSMA High:",   nrow(sig_df %>% filter(log2FoldChange > 0))),
               paste("Down PSMA High:", nrow(sig_df %>% filter(log2FoldChange < 0)))),
             file.path(cohort_dir, paste0(cohort_label, "_DESeq2_summary.txt")))
  
  run_go_enrichment(sig_df = sig_df, gene_map = gene_map,
                    cohort_dir = cohort_dir, prefix = cohort_label,
                    positive_label = "PSMA High", negative_label = "PSMA Low")
}

# -----------------------------------------------------------------------------
# 8. Helpers survie
# -----------------------------------------------------------------------------

run_group_survival_analysis <- function(analysis_df, group_var, group_levels,
                                        cohort_dir, prefix, title_text,
                                        normalize_plot_time = FALSE,
                                        xlab_text = "OS (days)") {
  keep    <- !is.na(analysis_df$OS_days) & analysis_df$OS_days > 0 &
    !is.na(analysis_df$death_event) & !is.na(analysis_df[[group_var]])
  surv_df <- analysis_df[keep, , drop = FALSE]
  write.csv(surv_df,
            file.path(cohort_dir, paste0(prefix, "_survival_input.csv")),
            row.names = FALSE)
  
  if (nrow(surv_df) < 4 || length(unique(surv_df[[group_var]])) < 2) {
    writeLines("Survival skipped: insufficient samples.",
               file.path(cohort_dir, paste0(prefix, "_survival_summary.txt")))
    return(invisible(NULL))
  }
  
  surv_df[[group_var]] <- factor(surv_df[[group_var]], levels = group_levels)
  sf  <- make_survival_formula(group_var)
  fit <- survfit(sf, data = surv_df)
  cox <- coxph(sf, data = surv_df)
  fit$call$formula <- sf; cox$call$formula <- sf
  cs  <- summary(cox)
  
  capture.output(cs, file = file.path(cohort_dir, paste0(prefix, "_cox_summary.txt")))
  write.csv(data.frame(
    variable    = rownames(cs$coefficients),
    HR          = cs$coefficients[, "exp(coef)"],
    lower_95_CI = cs$conf.int[, "lower .95"],
    upper_95_CI = cs$conf.int[, "upper .95"],
    p_value     = cs$coefficients[, "Pr(>|z|)"],
    row.names   = NULL),
    file.path(cohort_dir, paste0(prefix, "_cox_hr.csv")), row.names = FALSE)
  
  pal <- c("Initial"  = "#4393C3", "Recurrent" = "#D6604D",
           "PSMA Low" = "#4393C3", "PSMA High" = "#D6604D")

  plot_df    <- surv_df
  plot_fit   <- fit
  plot_title <- title_text
  if (normalize_plot_time) {
    origin_df <- surv_df %>%
      group_by(.data[[group_var]]) %>%
      summarise(OS_origin_days = min(OS_days, na.rm = TRUE),
                n_samples = dplyr::n(), .groups = "drop")
    plot_df <- surv_df %>%
      left_join(origin_df, by = group_var) %>%
      mutate(OS_days_original = OS_days,
             OS_days = pmax(OS_days - OS_origin_days, 0))
    write.csv(origin_df,
              file.path(cohort_dir, paste0(prefix, "_survival_time0_offsets.csv")),
              row.names = FALSE)
    write.csv(plot_df,
              file.path(cohort_dir, paste0(prefix, "_survival_input_time0_normalized.csv")),
              row.names = FALSE)
    plot_fit <- survfit(sf, data = plot_df)
    plot_fit$call$formula <- sf
    plot_title <- paste0(title_text, " (time normalized)")
  }

  p_km <- ggsurvplot(plot_fit, data = plot_df, risk.table = TRUE, pval = TRUE,
                     conf.int = FALSE, palette = unname(pal[group_levels]),
                     title = plot_title, xlab = xlab_text,
                     ylab = "OS probability", axes.offset = FALSE)
  if (normalize_plot_time) {
    p_km$plot <- p_km$plot +
      labs(subtitle = paste0("Time zero rebased within each ", group_var,
                             " group to its earliest observed OS day"))
  }
  save_plot(file.path(cohort_dir, paste0(prefix, "_KM.png")),
            p_km$plot, width = 6.8, height = 5.2)
}

run_survival_analysis <- function(df, cohort_dir, cohort_label) {
  normalize_recurrent <- identical(cohort_label, "Recurrent")
  run_group_survival_analysis(df, "PSMA_group", c("PSMA Low", "PSMA High"),
                              cohort_dir, cohort_label,
                              paste0(cohort_label, " OS: PSMA High vs Low"),
                              normalize_plot_time = normalize_recurrent,
                              xlab_text = if (normalize_recurrent)
                                "OS from recurrent time zero (days)"
                              else "OS (days)")
}

# -----------------------------------------------------------------------------
# 9. run_angiogenesis_volcano
# -----------------------------------------------------------------------------

run_angiogenesis_volcano <- function(res_shrunk_df, angio_ids, angio_symbols,
                                     cohort_dir, prefix) {
  if (is.null(res_shrunk_df) || nrow(res_shrunk_df) == 0) return(invisible(NULL))
  
  vdf <- res_shrunk_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      label    = dplyr::if_else(is.na(SYMBOL) | SYMBOL == "", gene_id, SYMBOL),
      nlp      = -log10(pmax(padj, 1e-300)),
      is_angio = gene_id %in% angio_ids,
      cat      = dplyr::case_when(
        is_angio & padj < 0.05 & log2FoldChange >  0 ~ "Angio - Up in Recurrent",
        is_angio & padj < 0.05 & log2FoldChange <= 0 ~ "Angio - Up in Initial",
        is_angio                                      ~ "Angio - NS",
        TRUE                                          ~ "Other"))
  
  n_det <- sum(vdf$is_angio)
  n_sig <- sum(vdf$is_angio & vdf$padj < 0.05)
  
  if (n_det == 0) {
    writeLines("Angiogenesis volcano skipped: no genes detected.",
               file.path(cohort_dir, paste0(prefix, "_angio_volcano_summary.txt")))
    return(invisible(NULL))
  }
  
  write.csv(vdf %>% filter(is_angio) %>%
              select(gene_id, SYMBOL, log2FoldChange, pvalue, padj, cat),
            file.path(cohort_dir, paste0(prefix, "_angiogenesis_genes_DESeq2.csv")),
            row.names = FALSE)
  
  col_vals   <- c("Angio - Up in Recurrent" = "#B2182B",
                  "Angio - Up in Initial"   = "#2166AC",
                  "Angio - NS"              = "#F4A582",
                  "Other"                   = "grey80")
  size_vals  <- c("Angio - Up in Recurrent" = 3, "Angio - Up in Initial" = 3,
                  "Angio - NS" = 2.5, "Other" = 1.0)
  alpha_vals <- c("Angio - Up in Recurrent" = 1, "Angio - Up in Initial" = 1,
                  "Angio - NS" = 0.8, "Other" = 0.3)
  
  p_a <- ggplot(vdf, aes(log2FoldChange, nlp,
                         color = cat, size = cat, alpha = cat)) +
    geom_point() +
    geom_label_repel(data = vdf %>% filter(is_angio),
                     aes(label = label), size = 2.8, max.overlaps = 60,
                     show.legend = FALSE, box.padding = 0.35,
                     segment.color = "grey50", segment.size = 0.3) +
    scale_color_manual(values = col_vals, name = "") +
    scale_size_manual(values = size_vals, guide = "none") +
    scale_alpha_manual(values = alpha_vals, guide = "none") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               linewidth = 0.4, color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.4, color = "grey40") +
    theme_classic(base_size = 13) +
    theme(legend.position = "top",
          plot.subtitle   = element_text(size = 9, color = "grey40")) +
    labs(title    = paste0(prefix, " — Angiogenesis gene signature"),
         subtitle = sprintf(
           "HALLMARK_ANGIOGENESIS: %d/%d genes detected, %d significant (FDR<0.05)",
           n_det, length(angio_ids), n_sig),
         x = "log2FC (Recurrent vs Initial; apeglm)", y = "-log10 adj. p")
  save_plot(file.path(cohort_dir, paste0(prefix, "_angiogenesis_volcano.png")),
            p_a, width = 8.5, height = 6.5)
  
  writeLines(c(paste("Prefix:", prefix),
               paste("Gene set size:", length(angio_ids)),
               paste("Detected:", n_det),
               paste("Sig FDR<0.05:", n_sig),
               paste("Up Recurrent:", sum(vdf$is_angio & vdf$padj < 0.05 & vdf$log2FoldChange > 0)),
               paste("Up Initial:",   sum(vdf$is_angio & vdf$padj < 0.05 & vdf$log2FoldChange < 0))),
             file.path(cohort_dir, paste0(prefix, "_angio_volcano_summary.txt")))
}

# -----------------------------------------------------------------------------
# 10. prepare_clinical_for_psm
# Actual columns in the clinical file: bcr_patient_barcode, gender (uppercase),
# age_at_initial_pathologic_diagnosis, pharmaceutical_tx_adjuvant,
# radiation_treatment_adjuvant
# -----------------------------------------------------------------------------

prepare_clinical_for_psm <- function(clinical_df) {
  
  psm_df <- clinical_df %>%
    transmute(
      case_id = bcr_patient_barcode,
      
      age_years = safe_numeric(age_at_initial_pathologic_diagnosis),
      
      gender_binary = dplyr::case_when(
        tolower(gender) == "male"   ~ 1L,
        tolower(gender) == "female" ~ 0L,
        TRUE ~ NA_integer_),
      
      pharma_treatment = dplyr::case_when(
        tolower(clean_tcga_val(pharmaceutical_tx_adjuvant)) == "yes" ~ 1L,
        tolower(clean_tcga_val(pharmaceutical_tx_adjuvant)) == "no"  ~ 0L,
        TRUE ~ NA_integer_),
      
      radiation_treatment = dplyr::case_when(
        tolower(clean_tcga_val(radiation_treatment_adjuvant)) == "yes" ~ 1L,
        tolower(clean_tcga_val(radiation_treatment_adjuvant)) == "no"  ~ 0L,
        TRUE ~ NA_integer_)
    )
  
  message("  tumor_grade absent from clinical file -> excluded from PSM.")
  
  na_pct <- psm_df %>%
    summarise(across(where(is.numeric), ~ round(mean(is.na(.)) * 100, 1)))
  message("  % NA per PSM covariate:")
  print(na_pct)
  
  psm_df
}

# -----------------------------------------------------------------------------
# 11. run_propensity_matching
# -----------------------------------------------------------------------------

run_propensity_matching <- function(counts_meta, clinical_psm_df, psm_dir,
                                    ratio = 3L) {
  ensure_dir(psm_dir)
  
  recurrent_patients    <- counts_meta %>%
    filter(sample_type_code == "02") %>% pull(patient_id) %>% unique()
  initial_only_patients <- counts_meta %>%
    filter(sample_type_code == "01", !patient_id %in% recurrent_patients) %>%
    pull(patient_id) %>% unique()
  
  psm_input <- clinical_psm_df %>%
    filter(case_id %in% c(recurrent_patients, initial_only_patients)) %>%
    mutate(is_recurrent = as.integer(case_id %in% recurrent_patients)) %>%
    filter(!is.na(age_years), !is.na(gender_binary))
  
  n_rec <- sum(psm_input$is_recurrent == 1)
  n_ini <- sum(psm_input$is_recurrent == 0)
  message(sprintf("  PSM input: %d recurrent, %d initial-only patients.", n_rec, n_ini))
  
  if (n_rec < 2 || n_ini < ratio) {
    writeLines(sprintf("PSM skipped: recurrent=%d, initial=%d.", n_rec, n_ini),
               file.path(psm_dir, "PSM_skipped.txt"))
    return(NULL)
  }
  
  # Covariates: keep only those with sufficient non-NA values
  covariates <- "age_years + gender_binary"
  if (mean(!is.na(psm_input$pharma_treatment)) > 0.5)
    covariates <- paste(covariates, "+ pharma_treatment")
  if (mean(!is.na(psm_input$radiation_treatment)) > 0.5)
    covariates <- paste(covariates, "+ radiation_treatment")
  
  psm_formula <- as.formula(paste("is_recurrent ~", covariates))
  message(sprintf("  PSM formula: %s", deparse(psm_formula)))
  
  set.seed(42)
  m_out <- tryCatch(
    matchit(formula = psm_formula, data = psm_input,
            method = "nearest", distance = "glm",
            ratio = ratio, replace = FALSE),
    error = function(e) { message("  MatchIt error: ", e$message); NULL })
  
  if (is.null(m_out)) {
    writeLines("PSM failed.", file.path(psm_dir, "PSM_failed.txt"))
    return(NULL)
  }
  
  matched <- match.data(m_out)
  ini_pts  <- matched %>% filter(is_recurrent == 0) %>% pull(case_id)
  rec_pts  <- matched %>% filter(is_recurrent == 1) %>% pull(case_id)
  ini_bc   <- counts_meta %>%
    filter(patient_id %in% ini_pts, sample_type_code == "01") %>%
    pull(sample_barcode)
  rec_bc   <- counts_meta %>%
    filter(patient_id %in% rec_pts, sample_type_code == "02") %>%
    pull(sample_barcode)
  
  write.csv(matched,
            file.path(psm_dir, "PSM_matched_patients.csv"), row.names = FALSE)
  write.csv(data.frame(sample_barcode = ini_bc, group = "Initial"),
            file.path(psm_dir, "PSM_matched_initial_barcodes.csv"), row.names = FALSE)
  write.csv(data.frame(sample_barcode = rec_bc, group = "Recurrent"),
            file.path(psm_dir, "PSM_matched_recurrent_barcodes.csv"), row.names = FALSE)
  writeLines(capture.output(summary(m_out, interactions = FALSE)),
             file.path(psm_dir, "PSM_balance_summary.txt"))
  
  tryCatch({
    # Love plot (ggplot) - much more readable than MatchIt's base R plot
    psm_sum   <- summary(m_out, interactions = FALSE)
    sum_df    <- as.data.frame(psm_sum$sum.all)
    sum_df$covariate <- rownames(sum_df)
    
    # Retrieve SMD before and after matching
    smd_all     <- abs(sum_df[["Std. Mean Diff."]])
    sum_matched <- as.data.frame(psm_sum$sum.matched)
    sum_matched$covariate <- rownames(sum_matched)
    smd_matched <- abs(sum_matched[["Std. Mean Diff."]])
    
    # Human-readable labels for covariates
    readable_names <- c(
      distance         = "Propensity score",
      age_years        = "Age (years)",
      gender_binary    = "Gender (Male=1)",
      pharma_treatment = "Chemo adjuvant",
      radiation_treatment = "Radiotherapy"
    )
    labels <- dplyr::coalesce(readable_names[sum_df$covariate], sum_df$covariate)
    
    love_df <- data.frame(
      covariate = factor(labels, levels = rev(labels)),
      Before    = smd_all,
      After     = smd_matched
    ) %>%
      tidyr::pivot_longer(cols = c("Before", "After"),
                          names_to = "Timing", values_to = "SMD") %>%
      mutate(Timing = factor(Timing, levels = c("Before", "After")))
    
    p_love <- ggplot(love_df, aes(x = SMD, y = covariate,
                                  color = Timing, shape = Timing)) +
      geom_vline(xintercept = 0.10, linetype = "dashed",
                 color = "grey40", linewidth = 0.6) +
      geom_vline(xintercept = 0.05, linetype = "dotted",
                 color = "grey60", linewidth = 0.5) +
      geom_point(size = 3.5, alpha = 0.9) +
      scale_color_manual(values = c("Before" = "#9ECAE1", "After" = "#2171B5"),
                         name = "Sample set") +
      scale_shape_manual(values = c("Before" = 1, "After" = 16),
                         name = "Sample set") +
      scale_x_continuous(limits = c(0, max(love_df$SMD, na.rm = TRUE) * 1.15),
                         breaks = c(0, 0.05, 0.10, 0.20, 0.30, 0.40)) +
      theme_classic(base_size = 13) +
      theme(
        legend.position  = "bottom",
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        axis.text.y      = element_text(size = 11)
      ) +
      labs(
        title    = sprintf("PSM Balance — %d:1 nearest-neighbour matching", ratio),
        subtitle = paste0("Dashed line = 0.10, Dotted = 0.05 (acceptable balance thresholds)\n",
                          "Lower SMD after matching = better covariate balance"),
        x  = "Absolute Standardised Mean Difference (SMD)",
        y  = NULL
      )
    
    ggsave(file.path(psm_dir, "PSM_balance_plot.png"),
           p_love, width = 7, height = 4, dpi = 300)
    
    # Text summary of covariate balance
    balance_ok <- love_df %>%
      filter(Timing == "After") %>%
      mutate(balanced = SMD < 0.10) %>%
      group_by(covariate) %>%
      summarise(SMD_after = round(SMD, 3), balanced = balanced, .groups = "drop")
    message("  PSM balance (SMD after matching):")
    print(balance_ok)
    
  }, error = function(e) {
    message("  PSM love plot failed: ", e$message)
  })
  
  writeLines(c(paste("Method: nearest neighbor,", ratio, ":1"),
               paste("Covariates:", covariates),
               paste("Matched recurrent patients:", length(rec_pts)),
               paste("Matched initial patients:", length(ini_pts)),
               paste("Matched initial barcodes:", length(ini_bc)),
               paste("Matched recurrent barcodes:", length(rec_bc))),
             file.path(psm_dir, "PSM_summary.txt"))
  
  message(sprintf("  PSM done: %d initial matches for %d recurrent.",
                  length(ini_pts), length(rec_pts)))
  list(initial_barcodes   = ini_bc,
       recurrent_barcodes = rec_bc,
       matched_patients   = matched)
}

# -----------------------------------------------------------------------------
# 12. analyse_cohort
# -----------------------------------------------------------------------------

analyse_cohort <- function(cohort_code, cohort_label,
                           counts_mat, counts_meta, fpkm_mat, fpkm_meta,
                           survival_df, hallmark_hypoxia_ids,
                           hallmark_hypoxia_symbols, gene_map_all, gene_map_fpkm) {
  
  cohort_dir <- file.path(output_root, cohort_label)
  ensure_dir(cohort_dir)
  
  cohort_samples <- counts_meta %>%
    filter(sample_type_code == cohort_code) %>% pull(sample_barcode)
  
  if (length(cohort_samples) == 0) {
    writeLines(paste0("No samples for cohort ", cohort_label, "."),
               file.path(cohort_dir, paste0(cohort_label, "_status.txt")))
    return(invisible(NULL))
  }
  
  counts_sub <- counts_mat[, cohort_samples, drop = FALSE]
  counts_sub <- counts_sub[rowSums(counts_sub) > 0, , drop = FALSE]
  
  dds_norm <- DESeqDataSetFromMatrix(countData = round(counts_sub),
                                     colData   = data.frame(row.names = cohort_samples),
                                     design    = ~ 1)
  dds_norm <- dds_norm[rowSums(counts(dds_norm)) > 10, ]
  dds_norm <- estimateSizeFactors(dds_norm)
  norm_counts <- counts(dds_norm, normalized = TRUE)
  vsd         <- vst(dds_norm, blind = TRUE)
  vsd_mat     <- assay(vsd)
  
  if (!folh1_id %in% rownames(vsd_mat))
    stop("FOLH1 (ENSG00000086205) absent after filtering.")
  
  cohort_meta <- counts_meta %>%
    filter(sample_barcode %in% colnames(vsd_mat)) %>%
    arrange(match(sample_barcode, colnames(vsd_mat))) %>%
    mutate(
      FOLH1_vst             = as.numeric(vsd_mat[folh1_id, sample_barcode]),
      FOLH1_norm_count      = as.numeric(norm_counts[folh1_id, sample_barcode]),
      FOLH1_log2_norm_count = log2(FOLH1_norm_count + 1),
      PSMA_group            = make_psma_groups(FOLH1_vst)) %>%
    left_join(survival_df, by = c("patient_id" = "case_id"))
  rownames(cohort_meta) <- cohort_meta$sample_barcode
  
  write.csv(cohort_meta,
            file.path(cohort_dir, paste0(cohort_label, "_sample_metadata.csv")),
            row.names = FALSE)
  write_gene_matrix(norm_counts, gene_map_all,
                    file.path(cohort_dir, paste0(cohort_label, "_normalized_counts.csv")))
  write_gene_matrix(vsd_mat, gene_map_all,
                    file.path(cohort_dir, paste0(cohort_label, "_VST_matrix.csv")))
  write.csv(cohort_meta %>%
              select(sample_barcode, patient_id, sample_type, PSMA_group,
                     FOLH1_vst, FOLH1_norm_count, FOLH1_log2_norm_count),
            file.path(cohort_dir, paste0(cohort_label, "_FOLH1_table.csv")),
            row.names = FALSE)
  
  p_f <- ggplot(cohort_meta,
                aes(PSMA_group, FOLH1_log2_norm_count, fill = PSMA_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, alpha = 0.55, size = 1.6) +
    scale_fill_manual(values = c("PSMA Low" = "#4393C3", "PSMA High" = "#D6604D")) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.margin = margin(t = 20, r = 10, b = 10, l = 10)) +
    labs(title = paste0(cohort_label, " FOLH1 by PSMA group"),
         x = "", y = "log2(DESeq2 norm. FOLH1 + 1)")
  p_f <- add_ttest_annotation(p_f, cohort_meta, "FOLH1_log2_norm_count",
                              "PSMA_group", "PSMA Low", "PSMA High")
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_FOLH1_boxplot.png")),
            p_f, width = 5.4, height = 4.6)
  
  gc <- cohort_meta %>% dplyr::count(PSMA_group) %>%
    mutate(p = n / sum(n), lbl = paste0(n, " (", round(100*p, 1), "%)"))
  p_b <- ggplot(gc, aes(PSMA_group, p, fill = PSMA_group)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = lbl), vjust = -0.4, size = 4) +
    scale_fill_manual(values = c("PSMA Low" = "#4393C3", "PSMA High" = "#D6604D")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, max(gc$p) * 1.2)) +
    theme_classic(base_size = 13) + theme(legend.position = "none") +
    labs(title = paste0(cohort_label, " PSMA groups"), x = "", y = "Proportion")
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_PSMA_proportions.png")),
            p_b, width = 5.2, height = 4.2)
  
  run_survival_analysis(cohort_meta, cohort_dir, cohort_label)
  
  run_hypoxia_analysis(vsd_mat, cohort_meta,
                       hallmark_hypoxia_ids, hallmark_hypoxia_symbols,
                       cohort_dir, cohort_label)
  
  if (!is.null(fpkm_mat) && !is.null(fpkm_meta)) {
    fpkm_bc     <- fpkm_meta %>%
      filter(sample_type_code == cohort_code) %>% pull(sample_barcode)
    only_fpkm   <- setdiff(fpkm_bc, cohort_meta$sample_barcode)
    only_counts <- setdiff(cohort_meta$sample_barcode, fpkm_bc)
    if (length(only_fpkm) > 0 || length(only_counts) > 0) {
      writeLines(c(paste("In FPKM only:", length(only_fpkm)),
                   if (length(only_fpkm) > 0) paste0("  ", only_fpkm),
                   paste("In counts only:", length(only_counts)),
                   if (length(only_counts) > 0) paste0("  ", only_counts)),
                 file.path(cohort_dir, paste0(cohort_label, "_barcode_mismatch.txt")))
    }
    fpkm_meta_sub <- fpkm_meta %>%
      filter(sample_type_code == cohort_code,
             sample_barcode %in% cohort_meta$sample_barcode) %>%
      left_join(cohort_meta %>%
                  select(sample_barcode, patient_id, sample_type,
                         PSMA_group, FOLH1_vst),
                by = c("sample_barcode", "patient_id", "sample_type"))
    run_fpkm_selected_gene_analysis(fpkm_mat, fpkm_meta_sub, gene_map_fpkm,
                                    cohort_dir, cohort_label)
  }
  
  dge_meta <- cohort_meta %>%
    select(sample_barcode, patient_id, sample_type,
           sample_type_code, sample_key, FOLH1_vst, PSMA_group)
  rownames(dge_meta) <- dge_meta$sample_barcode
  run_deseq_analysis(counts_sub[, dge_meta$sample_barcode, drop = FALSE],
                     dge_meta, gene_map_all, cohort_dir, cohort_label)
  
  writeLines(c(paste("Cohort:", cohort_label),
               paste("Samples:", nrow(cohort_meta)),
               paste("Patients:", dplyr::n_distinct(cohort_meta$patient_id)),
               paste("PSMA Low:", sum(cohort_meta$PSMA_group == "PSMA Low")),
               paste("PSMA High:", sum(cohort_meta$PSMA_group == "PSMA High")),
               paste("Median FOLH1 VST:",
                     round(median(cohort_meta$FOLH1_vst, na.rm = TRUE), 4))),
             file.path(cohort_dir, paste0(cohort_label, "_cohort_summary.txt")))
}

# -----------------------------------------------------------------------------
# 13. run_initial_vs_recurrent_analysis
# -----------------------------------------------------------------------------

run_initial_vs_recurrent_analysis <- function(
    counts_mat, counts_meta, fpkm_mat, fpkm_meta,
    survival_df, hallmark_hypoxia_ids, hallmark_hypoxia_symbols,
    angiogenesis_ids, angiogenesis_symbols,
    gene_map_all, gene_map_fpkm,
    psm_result = NULL) {
  
  cohort_dir <- file.path(output_root, "Initial_vs_Recurrent")
  ensure_dir(cohort_dir)
  
  if (!is.null(psm_result)) {
    ini_bc <- psm_result$initial_barcodes
    rec_bc <- psm_result$recurrent_barcodes
    writeLines(paste("Mode: PSM-matched initial (n =", length(ini_bc),
                     ") vs all recurrent (n =", length(rec_bc), ")"),
               file.path(cohort_dir, "sample_selection.txt"))
  } else {
    ini_bc <- counts_meta %>% filter(sample_type_code == "01") %>% pull(sample_barcode)
    rec_bc <- counts_meta %>% filter(sample_type_code == "02") %>% pull(sample_barcode)
    writeLines("Mode: all initial vs all recurrent (no PSM).",
               file.path(cohort_dir, "sample_selection.txt"))
  }
  
  comp_samples <- c(ini_bc, rec_bc)
  if (length(comp_samples) == 0) {
    writeLines("No samples.", file.path(cohort_dir, "status.txt"))
    return(invisible(NULL))
  }
  
  counts_sub <- counts_mat[, comp_samples, drop = FALSE]
  counts_sub <- counts_sub[rowSums(counts_sub) > 0, , drop = FALSE]
  
  dds_norm <- DESeqDataSetFromMatrix(countData = round(counts_sub),
                                     colData   = data.frame(row.names = comp_samples),
                                     design    = ~ 1)
  dds_norm <- dds_norm[rowSums(counts(dds_norm)) > 10, ]
  dds_norm <- estimateSizeFactors(dds_norm)
  norm_counts <- counts(dds_norm, normalized = TRUE)
  vsd         <- vst(dds_norm, blind = TRUE)
  vsd_mat     <- assay(vsd)
  
  comp_meta <- counts_meta %>%
    filter(sample_barcode %in% colnames(vsd_mat),
           sample_barcode %in% comp_samples) %>%
    arrange(match(sample_barcode, colnames(vsd_mat))) %>%
    mutate(
      Comparison_group      = factor(sample_type, levels = c("Initial", "Recurrent")),
      FOLH1_vst             = as.numeric(vsd_mat[folh1_id, sample_barcode]),
      FOLH1_norm_count      = as.numeric(norm_counts[folh1_id, sample_barcode]),
      FOLH1_log2_norm_count = log2(FOLH1_norm_count + 1)) %>%
    left_join(survival_df, by = c("patient_id" = "case_id"))
  rownames(comp_meta) <- comp_meta$sample_barcode
  
  ini_sub <- comp_meta %>% filter(Comparison_group == "Initial")
  if (nrow(ini_sub) >= 2 && length(unique(ini_sub$FOLH1_vst)) >= 2) {
    ini_groups <- make_psma_groups(ini_sub$FOLH1_vst)
    comp_meta$PSMA_group_matched <- NA_character_
    comp_meta$PSMA_group_matched[comp_meta$Comparison_group == "Initial"] <-
      as.character(ini_groups)
    message(sprintf("  PSMA H/L redefini dans %d initial matches.", nrow(ini_sub)))
  }
  
  write.csv(comp_meta,
            file.path(cohort_dir, "Initial_vs_Recurrent_metadata.csv"),
            row.names = FALSE)
  
  title_folh1 <- if (!is.null(psm_result)) "FOLH1: Matched Initial vs Recurrent"
  else "FOLH1: Initial vs Recurrent"
  p_f <- ggplot(comp_meta,
                aes(Comparison_group, FOLH1_log2_norm_count, fill = Comparison_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, alpha = 0.55, size = 1.6) +
    scale_fill_manual(values = c("Initial" = "#4393C3", "Recurrent" = "#D6604D")) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.margin = margin(t = 20, r = 10, b = 10, l = 10)) +
    labs(title = title_folh1, x = "", y = "log2(DESeq2 norm. FOLH1 + 1)")
  p_f <- add_ttest_annotation(p_f, comp_meta, "FOLH1_log2_norm_count",
                              "Comparison_group", "Initial", "Recurrent")
  save_plot(file.path(cohort_dir, "Initial_vs_Recurrent_FOLH1_boxplot.png"),
            p_f, width = 5.4, height = 4.6)
  
  gc <- comp_meta %>% dplyr::count(Comparison_group) %>%
    mutate(p = n / sum(n), lbl = paste0(n, " (", round(100*p, 1), "%)"))
  p_b <- ggplot(gc, aes(Comparison_group, p, fill = Comparison_group)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = lbl), vjust = -0.4, size = 4) +
    scale_fill_manual(values = c("Initial" = "#4393C3", "Recurrent" = "#D6604D")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, max(gc$p) * 1.2)) +
    theme_classic(base_size = 13) + theme(legend.position = "none") +
    labs(title = "Sample proportions", x = "", y = "Proportion")
  save_plot(file.path(cohort_dir, "Initial_vs_Recurrent_proportions.png"),
            p_b, width = 5.2, height = 4.2)
  
  run_group_survival_analysis(comp_meta, "Comparison_group",
                              c("Initial", "Recurrent"), cohort_dir,
                              "Initial_vs_Recurrent", "OS: Initial vs Recurrent")
  
  run_hypoxia_analysis(vsd_mat, comp_meta,
                       hallmark_hypoxia_ids, hallmark_hypoxia_symbols,
                       cohort_dir, "Initial_vs_Recurrent",
                       group_var    = "Comparison_group",
                       group_levels = c("Initial", "Recurrent"),
                       palette      = c("Initial" = "#4393C3", "Recurrent" = "#D6604D"))
  
  if (!is.null(fpkm_mat)) {
    common_fpkm <- intersect(comp_meta$sample_barcode, rownames(fpkm_mat))
    if (length(common_fpkm) >= 2) {
      fpkm_sub <- fpkm_mat[common_fpkm, , drop = FALSE]
      meta_sub  <- comp_meta %>%
        filter(sample_barcode %in% common_fpkm) %>%
        arrange(match(sample_barcode, common_fpkm))
      long_df <- as.data.frame(fpkm_sub) %>%
        rownames_to_column("sample_barcode") %>%
        tidyr::pivot_longer(-sample_barcode, names_to = "gene_id", values_to = "FPKM") %>%
        left_join(meta_sub %>% select(sample_barcode, Comparison_group, FOLH1_vst),
                  by = "sample_barcode") %>%
        left_join(gene_map_fpkm, by = c("gene_id" = "ENSEMBL")) %>%
        mutate(gene_label = ifelse(is.na(SYMBOL) | SYMBOL == "", gene_id, SYMBOL),
               log2_FPKM  = log2(FPKM + 1))
      stats_df <- long_df %>%
        group_by(gene_id, gene_label) %>%
        summarise(
          median_ini = median(FPKM[Comparison_group == "Initial"],   na.rm = TRUE),
          median_rec = median(FPKM[Comparison_group == "Recurrent"], na.rm = TRUE),
          log2FC     = log2(median_rec + 1) - log2(median_ini + 1),
          p_value    = { wt <- safe_wilcox_groups(FPKM, Comparison_group,
                                                  "Initial", "Recurrent")
          if (is.null(wt)) NA_real_ else wt$p.value },
          .groups = "drop") %>%
        mutate(padj = p.adjust(p_value, method = "BH")) %>%
        arrange(padj, desc(abs(log2FC)))
      write.csv(stats_df,
                file.path(cohort_dir, "Initial_vs_Recurrent_fpkm_stats.csv"),
                row.names = FALSE)
      heat_mat <- t(log2(fpkm_sub + 1))
      heat_gene_ids <- rownames(heat_mat)
      heat_labels <- gene_map_fpkm$SYMBOL[match(heat_gene_ids, gene_map_fpkm$ENSEMBL)]
      rownames(heat_mat) <- make.unique(ifelse(is.na(heat_labels) | heat_labels == "",
                                               heat_gene_ids, heat_labels))
      heat_z <- row_zscore(heat_mat)
      ordered_s <- cluster_order_from_z(heat_z, "cols")
      ordered_g <- cluster_order_from_z(heat_z[, ordered_s, drop = FALSE], "rows")
      hmap_df <- z_matrix_to_long(heat_z, sample_order = ordered_s,
                                  gene_order = ordered_g)
      p_h <- ggplot(hmap_df, aes(sample_barcode, gene_label, fill = z)) +
        geom_tile() +
        z_heatmap_fill() +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              panel.grid  = element_blank()) +
        labs(title = "Selected hypoxic genes (FPKM): Initial vs Recurrent",
             subtitle = "Rows and samples clustered by row z-score",
             x = "Samples", y = "Gene", fill = "z-score")
      save_plot(file.path(cohort_dir, "Initial_vs_Recurrent_fpkm_heatmap.png"),
                p_h, width = 8, height = 4.8)
    }
  }
  
  grp_n <- table(comp_meta$Comparison_group)
  if (length(grp_n) == 2 && all(grp_n >= 2)) {
    dds <- DESeqDataSetFromMatrix(
      countData = round(counts_sub[, comp_meta$sample_barcode, drop = FALSE]),
      colData   = comp_meta,
      design    = ~ Comparison_group)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- DESeq(dds)
    res        <- results(dds,
                          contrast = c("Comparison_group", "Recurrent", "Initial"),
                          alpha    = 0.05)
    res_shrunk <- lfcShrink(dds,
                            coef = get_coef_name(dds, "^Comparison_group_"),
                            type = "apeglm")
    ann <- function(r) as.data.frame(r) %>%
      rownames_to_column("gene_id") %>%
      left_join(gene_map_all, by = c("gene_id" = "ENSEMBL")) %>%
      arrange(padj)
    res_df        <- ann(res)
    res_shrunk_df <- ann(res_shrunk)
    sig_df        <- res_shrunk_df %>% filter(!is.na(padj), padj < 0.05)
    
    write.csv(res_df,
              file.path(cohort_dir, "Initial_vs_Recurrent_DESeq2_all.csv"),
              row.names = FALSE)
    write.csv(res_shrunk_df,
              file.path(cohort_dir, "Initial_vs_Recurrent_DESeq2_LFCshrink.csv"),
              row.names = FALSE)
    write.csv(sig_df,
              file.path(cohort_dir, "Initial_vs_Recurrent_DESeq2_sig.csv"),
              row.names = FALSE)
    
    make_volcano(
      df         = res_shrunk_df,
      up_label   = "Up in Recurrent",
      down_label = "Up in Initial",
      title      = "Volcano: Recurrent vs Initial",
      x_label    = "log2FC (Recurrent vs Initial; apeglm)",
      out_file   = file.path(cohort_dir, "Initial_vs_Recurrent_DESeq2_volcano.png")
    )
    
    run_angiogenesis_volcano(res_shrunk_df, angiogenesis_ids, angiogenesis_symbols,
                             cohort_dir, "Initial_vs_Recurrent")
    
    run_go_enrichment(sig_df = sig_df, gene_map = gene_map_all,
                      cohort_dir = cohort_dir, prefix = "Initial_vs_Recurrent",
                      positive_label = "Recurrent", negative_label = "Initial")
  } else {
    writeLines("DESeq2 skipped: need >= 2 samples per group.",
               file.path(cohort_dir, "Initial_vs_Recurrent_DESeq2_summary.txt"))
  }
}

# -----------------------------------------------------------------------------
# 14. run_paired_deseq_analysis
# Analyses the patients for whom both initial AND recurrent samples are available.
# DESeq2 design: ~ patient_id + sample_type  (removes inter-patient variability)
# -----------------------------------------------------------------------------

run_paired_deseq_analysis <- function(counts_mat, counts_meta,
                                      hallmark_hypoxia_ids, hallmark_hypoxia_symbols,
                                      angiogenesis_ids, angiogenesis_symbols,
                                      gene_map_all) {
  cohort_dir <- file.path(output_root, "Paired_Analysis")
  ensure_dir(cohort_dir)
  
  # ── Identify paired patients ────────────────────────────────────────────────
  initial_pts   <- counts_meta %>% filter(sample_type_code == "01") %>%
    pull(patient_id) %>% unique()
  recurrent_pts <- counts_meta %>% filter(sample_type_code == "02") %>%
    pull(patient_id) %>% unique()
  paired_pts    <- intersect(initial_pts, recurrent_pts)
  
  writeLines(c(
    paste("Paired patients found:", length(paired_pts)),
    if (length(paired_pts) > 0) paste0("  IDs: ", paste(paired_pts, collapse = ", "))
  ), file.path(cohort_dir, "paired_patients_list.txt"))
  
  if (length(paired_pts) < 2) {
    writeLines("Paired analysis skipped: need >= 2 paired patients.",
               file.path(cohort_dir, "paired_status.txt"))
    message("  Paired analysis skipped: only ", length(paired_pts), " paired patient(s).")
    return(invisible(NULL))
  }
  
  # One sample per patient per type (first occurrence if duplicates)
  paired_meta <- counts_meta %>%
    filter(patient_id %in% paired_pts, sample_type_code %in% c("01", "02")) %>%
    arrange(patient_id, sample_type_code) %>%
    group_by(patient_id, sample_type_code) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  # Keep only complete pairs
  complete_pts <- paired_meta %>% count(patient_id) %>%
    filter(n == 2) %>% pull(patient_id)
  paired_meta  <- paired_meta %>%
    filter(patient_id %in% complete_pts) %>%
    mutate(patient_id_f = factor(patient_id),
           sample_type  = factor(sample_type, levels = c("Initial", "Recurrent")))
  
  n_pairs <- length(complete_pts)
  message(sprintf("  Paired analysis: %d complete patient pairs (%d samples).",
                  n_pairs, nrow(paired_meta)))
  
  # ── VST normalisation (intercept-only, blind) ───────────────────────────────
  counts_sub <- counts_mat[, paired_meta$sample_barcode, drop = FALSE]
  counts_sub <- counts_sub[rowSums(counts_sub) > 0, , drop = FALSE]
  
  dds_norm <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData   = data.frame(row.names = paired_meta$sample_barcode),
    design    = ~ 1)
  dds_norm <- dds_norm[rowSums(counts(dds_norm)) > 10, ]
  dds_norm <- estimateSizeFactors(dds_norm)
  vsd      <- vst(dds_norm, blind = TRUE)
  vsd_mat  <- assay(vsd)
  
  paired_meta <- paired_meta %>%
    mutate(FOLH1_vst = as.numeric(vsd_mat[folh1_id, sample_barcode]))
  write.csv(paired_meta,
            file.path(cohort_dir, "paired_sample_metadata.csv"),
            row.names = FALSE)
  
  # ── DESeq2 paired design ────────────────────────────────────────────────────
  paired_folh1_wide <- paired_meta %>%
    select(patient_id, sample_type, FOLH1_vst) %>%
    tidyr::pivot_wider(names_from = sample_type, values_from = FOLH1_vst) %>%
    filter(!is.na(Initial), !is.na(Recurrent))
  paired_folh1_p <- if (nrow(paired_folh1_wide) >= 2) {
    tryCatch(t.test(paired_folh1_wide$Recurrent,
                    paired_folh1_wide$Initial,
                    paired = TRUE)$p.value,
             error = function(e) NA_real_)
  } else {
    NA_real_
  }
  write.csv(data.frame(
    n_pairs = nrow(paired_folh1_wide),
    mean_initial = mean(paired_folh1_wide$Initial, na.rm = TRUE),
    mean_recurrent = mean(paired_folh1_wide$Recurrent, na.rm = TRUE),
    median_initial = median(paired_folh1_wide$Initial, na.rm = TRUE),
    median_recurrent = median(paired_folh1_wide$Recurrent, na.rm = TRUE),
    paired_ttest_p = paired_folh1_p),
    file.path(cohort_dir, "Paired_FOLH1_paired_ttest.csv"),
    row.names = FALSE)
  
  folh1_y_max <- max(paired_meta$FOLH1_vst, na.rm = TRUE)
  folh1_y_min <- min(paired_meta$FOLH1_vst, na.rm = TRUE)
  folh1_y_ann <- folh1_y_max + (folh1_y_max - folh1_y_min) * 0.10
  p_folh1 <- ggplot(paired_meta,
                    aes(x = sample_type, y = FOLH1_vst, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.55) +
    geom_line(aes(group = patient_id), color = "grey45",
              linewidth = 0.35, alpha = 0.65) +
    geom_point(size = 2.2, alpha = 0.85, color = "black", shape = 21) +
    scale_fill_manual(values = c("Initial" = "#4393C3", "Recurrent" = "#D6604D")) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.margin = margin(t = 20, r = 10, b = 10, l = 10)) +
    annotate("text", x = 1.5, y = folh1_y_ann,
             label = paste0("paired ", format_pval(paired_folh1_p)),
             size = 4, fontface = "italic") +
    labs(title = sprintf("Paired FOLH1 expression (n=%d patients)",
                         nrow(paired_folh1_wide)),
         x = "", y = "FOLH1 VST expression")
  save_plot(file.path(cohort_dir, "Paired_FOLH1_expression_boxplot.png"),
            p_folh1, width = 5.4, height = 4.6)
  
  meta_deseq <- paired_meta %>%
    select(sample_barcode, patient_id_f, sample_type) %>%
    as.data.frame()
  rownames(meta_deseq) <- meta_deseq$sample_barcode
  
  counts_deseq <- round(counts_sub[, meta_deseq$sample_barcode, drop = FALSE])
  counts_deseq <- counts_deseq[rowSums(counts_deseq) >= 10, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_deseq,
    colData   = meta_deseq,
    design    = ~ patient_id_f + sample_type)
  dds <- DESeq(dds)
  
  res        <- results(dds,
                        contrast = c("sample_type", "Recurrent", "Initial"),
                        alpha    = 0.05)
  coef_name  <- get_coef_name(dds, "sample_type")
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  
  ann <- function(r) as.data.frame(r) %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_map_all, by = c("gene_id" = "ENSEMBL")) %>%
    arrange(padj)
  
  res_df        <- ann(res)
  res_shrunk_df <- ann(res_shrunk)
  sig_df        <- res_shrunk_df %>% filter(!is.na(padj), padj < 0.05)
  
  write.csv(res_df,
            file.path(cohort_dir, "Paired_DESeq2_all_genes.csv"),
            row.names = FALSE)
  write.csv(res_shrunk_df,
            file.path(cohort_dir, "Paired_DESeq2_LFCshrink.csv"),
            row.names = FALSE)
  write.csv(sig_df,
            file.path(cohort_dir, "Paired_DESeq2_sig_FDR0.05.csv"),
            row.names = FALSE)
  
  writeLines(c(
    paste("Design: ~ patient_id + sample_type (paired)"),
    paste("Paired patients:", n_pairs),
    paste("Total samples:", nrow(paired_meta)),
    paste("Genes tested:", nrow(res_df)),
    paste("Sig FDR<0.05:", nrow(sig_df)),
    paste("Up in Recurrent:", nrow(sig_df %>% filter(log2FoldChange > 0))),
    paste("Up in Initial:",   nrow(sig_df %>% filter(log2FoldChange < 0)))
  ), file.path(cohort_dir, "Paired_DESeq2_summary.txt"))
  
  # ── Volcano ──────────────────────────────────────────────────────────────────
  make_volcano(
    df         = res_shrunk_df,
    up_label   = "Up in Recurrent",
    down_label = "Up in Initial",
    title      = sprintf("Paired (n=%d patients) — Recurrent vs Initial", n_pairs),
    subtitle   = "Design: ~ patient_id + sample_type",
    x_label    = "log2FC (Recurrent vs Initial; apeglm)",
    out_file   = file.path(cohort_dir, "Paired_DESeq2_volcano.png")
  )
  
  # ── Angiogenesis volcano ────────────────────────────────────────────────────
  run_angiogenesis_volcano(res_shrunk_df, angiogenesis_ids, angiogenesis_symbols,
                           cohort_dir, "Paired")
  
  # ── Hypoxia score ────────────────────────────────────────────────────────────
  run_hypoxia_analysis(
    vsd_mat      = vsd_mat,
    cohort_meta  = paired_meta %>%
      mutate(PSMA_group = as.character(sample_type)),
    hallmark_ids     = hallmark_hypoxia_ids,
    hallmark_symbols = hallmark_hypoxia_symbols,
    cohort_dir       = cohort_dir,
    cohort_label     = "Paired",
    group_var        = "PSMA_group",
    group_levels     = c("Initial", "Recurrent"),
    palette          = c("Initial" = "#4393C3", "Recurrent" = "#D6604D")
  )
  
  # ── GO enrichment ───────────────────────────────────────────────────────────
  run_go_enrichment(sig_df = sig_df, gene_map = gene_map_all,
                    cohort_dir = cohort_dir, prefix = "Paired",
                    positive_label = "Recurrent", negative_label = "Initial")
  
  message(sprintf("  Paired analysis complete. %d DEGs (FDR<0.05).", nrow(sig_df)))
  invisible(sig_df)
}

# -----------------------------------------------------------------------------
# 14b. Clean plot overrides integrated from the standalone heatmap/volcano scripts
# -----------------------------------------------------------------------------

clean_volcano_base_mean_min <- 10
clean_volcano_lfc_cutoff    <- 1
clean_volcano_fdr_cutoff    <- 0.05

format_n_clean <- function(x) {
  format(x, big.mark = ",", scientific = FALSE)
}

prepare_volcano_data_clean <- function(df,
                                       lfc_col = "log2FoldChange",
                                       padj_col = "padj",
                                       label_col = "SYMBOL",
                                       gene_id_col = "gene_id",
                                       base_mean_min = clean_volcano_base_mean_min) {
  if (!label_col %in% colnames(df)) df[[label_col]] <- NA_character_
  if (!gene_id_col %in% colnames(df)) df[[gene_id_col]] <- rownames(df)
  
  vdf <- df %>%
    filter(!is.na(.data[[lfc_col]]),
           !is.na(.data[[padj_col]]),
           is.finite(.data[[lfc_col]]))
  
  if ("baseMean" %in% colnames(vdf)) {
    vdf <- vdf %>% filter(!is.na(baseMean), baseMean >= base_mean_min)
  }
  
  if (nrow(vdf) == 0) {
    stop("No genes left for volcano plotting after filtering.")
  }
  
  positive_padj <- vdf[[padj_col]][vdf[[padj_col]] > 0 & is.finite(vdf[[padj_col]])]
  padj_floor <- if (length(positive_padj) > 0) {
    max(min(positive_padj, na.rm = TRUE) * 0.5, .Machine$double.xmin)
  } else {
    .Machine$double.xmin
  }
  
  vdf %>%
    mutate(
      label_text = dplyr::if_else(
        is.na(.data[[label_col]]) | .data[[label_col]] == "",
        as.character(.data[[gene_id_col]]),
        as.character(.data[[label_col]])
      ),
      padj_plot = pmax(.data[[padj_col]], padj_floor),
      neglog10  = -log10(padj_plot),
      direction = dplyr::case_when(
        .data[[padj_col]] < clean_volcano_fdr_cutoff &
          .data[[lfc_col]] >= clean_volcano_lfc_cutoff ~ "up",
        .data[[padj_col]] < clean_volcano_fdr_cutoff &
          .data[[lfc_col]] <= -clean_volcano_lfc_cutoff ~ "down",
        TRUE ~ "ns"
      )
    )
}

choose_volcano_limits_clean <- function(vdf, lfc_col = "log2FoldChange") {
  max_abs_lfc <- max(abs(vdf[[lfc_col]]), na.rm = TRUE)
  x_lim <- ceiling(max(2, max_abs_lfc) * 10) / 10
  x_lim <- min(x_lim, 15)
  
  y_lim <- max(vdf$neglog10, na.rm = TRUE)
  y_lim <- ceiling(max(2, y_lim) * 10) / 10
  
  list(x = c(-x_lim, x_lim), y = c(0, y_lim * 1.08))
}

select_clean_volcano_labels <- function(vdf,
                                        up_label,
                                        down_label,
                                        lfc_col = "log2FoldChange",
                                        padj_col = "padj",
                                        n_each_side = 7,
                                        force_symbols = character(0)) {
  sig_hits <- vdf %>%
    filter(direction != "ns") %>%
    group_by(direction) %>%
    arrange(.data[[padj_col]], desc(abs(.data[[lfc_col]])), .by_group = TRUE) %>%
    slice_head(n = n_each_side) %>%
    ungroup()
  
  forced <- vdf %>%
    filter(label_text %in% force_symbols)
  
  bind_rows(sig_hits, forced) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    mutate(class = factor(dplyr::case_when(
      direction == "up"   ~ up_label,
      direction == "down" ~ down_label,
      TRUE                ~ "Not significant"
    ), levels = c(up_label, down_label, "Not significant")))
}

# Clean replacement for the old volcano helper. It keeps the same signature so
# all existing pipeline calls automatically use the improved plotting style.
make_volcano <- function(df,
                         lfc_col       = "log2FoldChange",
                         padj_col      = "padj",
                         label_col     = "SYMBOL",
                         gene_id_col   = "gene_id",
                         lfc_threshold = clean_volcano_lfc_cutoff,
                         padj_threshold = clean_volcano_fdr_cutoff,
                         up_label      = "Up",
                         down_label    = "Down",
                         n_labels      = 7,
                         title         = "",
                         subtitle      = "",
                         x_label       = "Shrunken log2 fold-change",
                         out_file      = NULL,
                         width = 7.4, height = 5.4) {
  
  old_lfc_cutoff <- clean_volcano_lfc_cutoff
  old_fdr_cutoff <- clean_volcano_fdr_cutoff
  clean_volcano_lfc_cutoff <<- lfc_threshold
  clean_volcano_fdr_cutoff <<- padj_threshold
  on.exit({
    clean_volcano_lfc_cutoff <<- old_lfc_cutoff
    clean_volcano_fdr_cutoff <<- old_fdr_cutoff
  }, add = TRUE)
  
  vdf <- prepare_volcano_data_clean(df, lfc_col, padj_col, label_col, gene_id_col)
  limits <- choose_volcano_limits_clean(vdf, lfc_col)
  
  vdf <- vdf %>%
    mutate(
      class = factor(dplyr::case_when(
        direction == "up"   ~ up_label,
        direction == "down" ~ down_label,
        TRUE                ~ "Not significant"
      ), levels = c(up_label, down_label, "Not significant"))
    )
  
  n_up   <- sum(vdf$direction == "up", na.rm = TRUE)
  n_down <- sum(vdf$direction == "down", na.rm = TRUE)
  clean_subtitle <- paste0(
    if (nchar(subtitle) > 0) paste0(subtitle, " | ") else "",
    "plotted genes: ", format_n_clean(nrow(vdf)),
    " (baseMean >= ", clean_volcano_base_mean_min, "); ",
    "FDR<", padj_threshold, " and |log2FC|>=", lfc_threshold,
    ": ", n_up, " up, ", n_down, " down"
  )
  
  force_symbols <- if (any(vdf$label_text == "FOLH1", na.rm = TRUE)) "FOLH1" else character(0)
  labels <- select_clean_volcano_labels(
    vdf,
    up_label = up_label,
    down_label = down_label,
    lfc_col = lfc_col,
    padj_col = padj_col,
    n_each_side = n_labels,
    force_symbols = force_symbols
  )
  
  cols <- c("#C84630", "#2C7FB8", "grey78")
  names(cols) <- c(up_label, down_label, "Not significant")
  
  p <- ggplot() +
    geom_point(
      data = vdf %>% filter(class == "Not significant"),
      aes(x = .data[[lfc_col]], y = neglog10),
      color = "grey72", alpha = 0.28, size = 0.85
    ) +
    geom_point(
      data = vdf %>% filter(class != "Not significant"),
      aes(x = .data[[lfc_col]], y = neglog10, color = class),
      alpha = 0.78, size = 1.15
    ) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_hline(yintercept = -log10(padj_threshold),
               linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_text_repel(
      data = labels,
      aes(x = .data[[lfc_col]], y = neglog10, label = label_text, color = class),
      size = 3.0,
      min.segment.length = 0,
      segment.size = 0.25,
      segment.alpha = 0.55,
      box.padding = 0.35,
      point.padding = 0.2,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(values = cols, drop = FALSE) +
    coord_cartesian(xlim = limits$x, ylim = limits$y, clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 8.8, color = "grey35"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 9.5),
      axis.title = element_text(size = 12.5),
      axis.text = element_text(size = 10.5),
      plot.margin = margin(t = 12, r = 18, b = 10, l = 10)
    ) +
    labs(title = title, subtitle = clean_subtitle,
         x = x_label, y = "-log10 adjusted p-value")
  
  if (!is.null(out_file)) {
    save_plot(out_file, p, width = width, height = height)
    write.csv(vdf, sub("\\.png$", "_plot_data_clean.csv", out_file),
              row.names = FALSE)
  }
  
  invisible(p)
}

make_folh1_top_heatmap <- function(vsd_mat_row_folh1,
                                   gene_z_mat,
                                   sample_order,
                                   gene_order,
                                   psma_groups = NULL,
                                   split_label = "PSMA Low",
                                   title = "",
                                   subtitle = "",
                                   x_label = "Samples",
                                   out_file = NULL,
                                   width = 9,
                                   height = 7) {
  n_samples <- length(sample_order)
  
  fv <- as.numeric(vsd_mat_row_folh1[1, sample_order])
  fmin <- min(fv, na.rm = TRUE)
  fmax <- max(fv, na.rm = TRUE)
  fz <- if (abs(fmax - fmin) > 1e-9) {
    (fv - fmin) / (fmax - fmin) * 5 - 2.5
  } else {
    rep(0, n_samples)
  }
  
  folh1_df <- data.frame(
    sample_barcode = factor(sample_order, levels = sample_order),
    gene = "FOLH1\n(min-max)",
    z = fz
  )
  
  do_split <- FALSE
  n_left <- NA_integer_
  if (!is.null(psma_groups) && !is.null(split_label)) {
    groups_in_order <- psma_groups[sample_order]
    n_left <- sum(groups_in_order == split_label, na.rm = TRUE)
    do_split <- n_left > 0 && n_left < n_samples
  }
  
  vline_layer <- if (do_split) {
    geom_vline(xintercept = n_left + 0.5,
               color = "black", linewidth = 0.8, linetype = "dashed")
  } else {
    NULL
  }
  
  p_top <- ggplot(folh1_df, aes(sample_barcode, gene, fill = z)) +
    geom_tile() +
    z_heatmap_fill("z-score") +
    vline_layer +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 4, r = 5, b = 0, l = 5)
    )
  
  gene_mat_ord <- gene_z_mat[gene_order, sample_order, drop = FALSE]
  gene_long <- data.frame(gene_label = rownames(gene_mat_ord),
                          gene_mat_ord, check.names = FALSE) %>%
    tidyr::pivot_longer(-gene_label,
                        names_to = "sample_barcode",
                        values_to = "z") %>%
    mutate(sample_barcode = factor(sample_barcode, levels = sample_order),
           gene_label = factor(gene_label, levels = rev(gene_order)))
  
  p_bottom <- ggplot(gene_long, aes(sample_barcode, gene_label, fill = z)) +
    geom_tile() +
    z_heatmap_fill("z-score") +
    vline_layer +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 7.5, color = "grey40"),
      plot.margin = margin(t = 0, r = 5, b = 10, l = 5)
    ) +
    labs(title = title, subtitle = subtitle, x = x_label, y = "Gene")
  
  p_out <- p_top / p_bottom +
    patchwork::plot_layout(heights = c(2, length(gene_order)),
                           guides = "collect") &
    theme(legend.position = "right")
  
  if (!is.null(out_file)) {
    save_plot(out_file, p_out, width = width, height = height)
  }
  
  invisible(p_out)
}

make_signature_heatmap_clean <- function(vsd_mat,
                                         cohort_meta,
                                         signature_ids,
                                         signature_symbols,
                                         cohort_dir,
                                         cohort_label,
                                         output_name,
                                         signature_label,
                                         group_var = "PSMA_group",
                                         split_label = "PSMA Low",
                                         top_n = 25) {
  common_samples <- intersect(cohort_meta$sample_barcode, colnames(vsd_mat))
  detected_ids <- intersect(signature_ids, rownames(vsd_mat))
  
  if (length(common_samples) < 2 || length(detected_ids) == 0 ||
      !folh1_id %in% rownames(vsd_mat)) {
    writeLines(sprintf("Heatmap skipped: samples=%d, detected genes=%d.",
                       length(common_samples), length(detected_ids)),
               file.path(cohort_dir, paste0(output_name, "_skipped.txt")))
    return(invisible(NULL))
  }
  
  z_mat <- row_zscore(vsd_mat[detected_ids, common_samples, drop = FALSE])
  sym_map <- map_ensembl_to_annotations(detected_ids)
  sym_labels <- sym_map$SYMBOL[match(detected_ids, sym_map$ENSEMBL)]
  rownames(z_mat) <- make.unique(ifelse(is.na(sym_labels) | sym_labels == "",
                                        detected_ids, sym_labels))
  
  sample_order <- cohort_meta %>%
    filter(sample_barcode %in% colnames(z_mat), !is.na(FOLH1_vst)) %>%
    arrange(FOLH1_vst) %>%
    pull(sample_barcode)
  if (length(sample_order) < 2) {
    sample_order <- cluster_order_from_z(z_mat, "cols")
  }
  
  gene_var <- apply(z_mat, 1, var, na.rm = TRUE)
  selected_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(top_n, length(gene_var)))]
  top_z <- z_mat[selected_genes, sample_order, drop = FALSE]
  gene_order <- cluster_order_from_z(top_z, "rows")
  
  group_vector <- NULL
  if (!is.null(group_var) && group_var %in% colnames(cohort_meta)) {
    group_vector <- setNames(as.character(cohort_meta[[group_var]]),
                             cohort_meta$sample_barcode)
  }
  
  subtitle <- paste0(
    "Samples ordered by FOLH1 low-to-high; FOLH1 strip is min-max scaled | ",
    length(detected_ids), "/", length(signature_symbols), " signature genes detected"
  )
  if (!is.null(split_label)) {
    subtitle <- paste0(subtitle, " | dashed line = ", split_label, " boundary")
  }
  
  make_folh1_top_heatmap(
    vsd_mat_row_folh1 = vsd_mat[folh1_id, sample_order, drop = FALSE],
    gene_z_mat = top_z,
    sample_order = sample_order,
    gene_order = gene_order,
    psma_groups = group_vector,
    split_label = split_label,
    title = paste0(cohort_label, " - FOLH1 + top ", length(gene_order),
                   " ", signature_label, " genes"),
    subtitle = subtitle,
    x_label = paste0("Samples (n=", length(sample_order), ")"),
    out_file = file.path(cohort_dir, paste0(output_name, ".png")),
    width = 9,
    height = 7
  )
}

# Clean replacement for the hypoxia workflow: keeps CSV/boxplot/scatter outputs,
# but regenerates the heatmap with the cleaner FOLH1 strip layout.
run_hypoxia_analysis <- function(vsd_mat, cohort_meta, hallmark_ids,
                                 hallmark_symbols, cohort_dir, cohort_label,
                                 group_var    = "PSMA_group",
                                 group_levels = c("PSMA Low", "PSMA High"),
                                 palette      = c("PSMA Low"  = "#4393C3",
                                                  "PSMA High" = "#D6604D")) {
  
  hypoxia_ids <- intersect(hallmark_ids, rownames(vsd_mat))
  if (length(hypoxia_ids) == 0) {
    writeLines("Hypoxia analysis skipped: no Hallmark genes in VST matrix.",
               file.path(cohort_dir, paste0(cohort_label, "_hypoxia_summary.txt")))
    return(invisible(NULL))
  }
  
  z_mat <- row_zscore(vsd_mat[hypoxia_ids, , drop = FALSE])
  hypoxia_score <- colMeans(z_mat, na.rm = TRUE)
  cohort_meta$Hallmark_Hypoxia_score <- hypoxia_score[cohort_meta$sample_barcode]
  
  hypoxia_df <- cohort_meta %>%
    select(sample_barcode, patient_id, sample_type,
           all_of(group_var), FOLH1_vst, Hallmark_Hypoxia_score)
  write.csv(hypoxia_df,
            file.path(cohort_dir, paste0(cohort_label, "_hypoxia_scores.csv")),
            row.names = FALSE)
  
  wt <- safe_wilcox(as.formula(paste("Hallmark_Hypoxia_score ~", group_var)),
                    hypoxia_df)
  grp_a <- group_levels[1]
  grp_b <- group_levels[2]
  write.csv(data.frame(
    cohort = cohort_label,
    n_hypoxia_genes = length(hypoxia_ids),
    group_a = grp_a,
    median_a = median(hypoxia_df$Hallmark_Hypoxia_score[hypoxia_df[[group_var]] == grp_a],
                      na.rm = TRUE),
    group_b = grp_b,
    median_b = median(hypoxia_df$Hallmark_Hypoxia_score[hypoxia_df[[group_var]] == grp_b],
                      na.rm = TRUE),
    wilcox_p = if (is.null(wt)) NA_real_ else wt$p.value),
    file.path(cohort_dir, paste0(cohort_label, "_hypoxia_group_stats.csv")),
    row.names = FALSE)
  
  sym_map <- map_ensembl_to_annotations(hypoxia_ids)
  gene_list_df <- data.frame(
    ENSEMBL = hypoxia_ids,
    SYMBOL = sym_map$SYMBOL[match(hypoxia_ids, sym_map$ENSEMBL)])
  write.csv(gene_list_df,
            file.path(cohort_dir, paste0(cohort_label, "_hypoxia_gene_list.csv")),
            row.names = FALSE)
  
  gene_syms_measured <- na.omit(gene_list_df$SYMBOL)
  subtitle_text <- paste0(
    "HALLMARK_HYPOXIA - ", length(hypoxia_ids), " genes measured",
    " (of ", length(hallmark_symbols), " in set) | ",
    "Top: ", paste(head(sort(gene_syms_measured), 10), collapse = ", "),
    if (length(gene_syms_measured) > 10) ", ..." else "")
  
  hypoxia_df[[group_var]] <- factor(hypoxia_df[[group_var]], levels = group_levels)
  
  p_box <- ggplot(hypoxia_df,
                  aes(x = .data[[group_var]], y = Hallmark_Hypoxia_score,
                      fill = .data[[group_var]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.55) +
    scale_fill_manual(values = palette) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = 7, color = "grey40"),
          plot.margin = margin(t = 20, r = 10, b = 10, l = 10)) +
    labs(title = paste0(cohort_label, " - Hallmark Hypoxia score"),
         subtitle = subtitle_text,
         x = "", y = "Hypoxia score (mean z-score)")
  p_box <- add_ttest_annotation(p_box, hypoxia_df, "Hallmark_Hypoxia_score",
                                group_var, group_levels[1], group_levels[2])
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_hypoxia_boxplot.png")),
            p_box, width = 6.5, height = 5.5)
  
  p_scatter <- ggplot(hypoxia_df,
                      aes(x = FOLH1_vst, y = Hallmark_Hypoxia_score,
                          color = .data[[group_var]])) +
    geom_point(size = 2.1, alpha = 0.8) +
    scale_color_manual(values = palette) +
    theme_classic(base_size = 13) +
    labs(title = paste0(cohort_label, " FOLH1 vs hypoxia score"),
         x = "FOLH1 VST expression", y = "Hallmark hypoxia score", color = "")
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_FOLH1_vs_hypoxia_scatter.png")),
            p_scatter, width = 5.8, height = 4.5)
  
  split_label <- if (identical(group_var, "PSMA_group")) "PSMA Low" else NULL
  make_signature_heatmap_clean(
    vsd_mat = vsd_mat,
    cohort_meta = cohort_meta,
    signature_ids = hallmark_ids,
    signature_symbols = hallmark_symbols,
    cohort_dir = cohort_dir,
    cohort_label = cohort_label,
    output_name = paste0(cohort_label, "_hypoxia_top25_heatmap"),
    signature_label = "Hallmark Hypoxia",
    group_var = group_var,
    split_label = split_label
  )
}

run_fpkm_selected_gene_analysis <- function(fpkm_mat, cohort_meta, gene_map,
                                            cohort_dir, cohort_label) {
  if (is.null(fpkm_mat) || nrow(fpkm_mat) == 0) {
    writeLines("FPKM analysis skipped.",
               file.path(cohort_dir, paste0(cohort_label, "_fpkm_summary.txt")))
    return(invisible(NULL))
  }
  
  common <- intersect(cohort_meta$sample_barcode, rownames(fpkm_mat))
  if (length(common) < 2) {
    writeLines("FPKM analysis skipped: < 2 overlapping samples.",
               file.path(cohort_dir, paste0(cohort_label, "_fpkm_summary.txt")))
    return(invisible(NULL))
  }
  
  meta_sub <- cohort_meta %>%
    filter(sample_barcode %in% common) %>%
    arrange(match(sample_barcode, common))
  fpkm_sub <- fpkm_mat[meta_sub$sample_barcode, , drop = FALSE]
  
  write.csv(as.data.frame(fpkm_sub) %>%
              rownames_to_column("sample_barcode") %>%
              left_join(meta_sub, by = "sample_barcode"),
            file.path(cohort_dir, paste0(cohort_label, "_fpkm_matrix.csv")),
            row.names = FALSE)
  
  long_df <- as.data.frame(fpkm_sub) %>%
    rownames_to_column("sample_barcode") %>%
    tidyr::pivot_longer(-sample_barcode, names_to = "gene_id", values_to = "FPKM") %>%
    left_join(meta_sub, by = "sample_barcode") %>%
    left_join(gene_map, by = c("gene_id" = "ENSEMBL")) %>%
    mutate(gene_label = ifelse(is.na(SYMBOL) | SYMBOL == "", gene_id, SYMBOL),
           log2_FPKM = log2(FPKM + 1))
  
  stats_df <- long_df %>%
    group_by(gene_id, gene_label) %>%
    summarise(
      median_low = median(FPKM[PSMA_group == "PSMA Low"], na.rm = TRUE),
      median_high = median(FPKM[PSMA_group == "PSMA High"], na.rm = TRUE),
      log2FC = log2(median_high + 1) - log2(median_low + 1),
      p_value = {
        wt <- safe_wilcox_groups(FPKM, PSMA_group, "PSMA Low", "PSMA High")
        if (is.null(wt)) NA_real_ else wt$p.value
      },
      .groups = "drop") %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    arrange(padj, desc(abs(log2FC)))
  write.csv(stats_df,
            file.path(cohort_dir, paste0(cohort_label, "_fpkm_stats.csv")),
            row.names = FALSE)
  
  heat_mat <- t(log2(fpkm_sub + 1))
  heat_gene_ids <- rownames(heat_mat)
  heat_labels <- gene_map$SYMBOL[match(heat_gene_ids, gene_map$ENSEMBL)]
  rownames(heat_mat) <- make.unique(ifelse(is.na(heat_labels) | heat_labels == "",
                                           heat_gene_ids, heat_labels))
  heat_z <- row_zscore(heat_mat)
  
  ordered_s <- meta_sub %>%
    filter(sample_barcode %in% colnames(heat_z), !is.na(FOLH1_vst)) %>%
    arrange(FOLH1_vst) %>%
    pull(sample_barcode)
  if (length(ordered_s) < 2) {
    ordered_s <- cluster_order_from_z(heat_z, "cols")
  }
  ordered_g <- cluster_order_from_z(heat_z[, ordered_s, drop = FALSE], "rows")
  hmap_df <- z_matrix_to_long(heat_z, sample_order = ordered_s,
                              gene_order = ordered_g)
  
  p_heat <- ggplot(hmap_df, aes(sample_barcode, gene_label, fill = z)) +
    geom_tile() +
    z_heatmap_fill() +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          plot.subtitle = element_text(size = 8, color = "grey40")) +
    labs(title = paste0(cohort_label, " - selected genes (FPKM)"),
         subtitle = "Samples ordered by FOLH1 expression low-to-high; genes clustered by z-score",
         x = paste0("Samples (n=", length(ordered_s), ")"),
         y = "Gene", fill = "z-score")
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_fpkm_heatmap.png")),
            p_heat, width = 8, height = 4.8)
}

run_angiogenesis_psma_analysis <- function(vsd_mat, cohort_meta, angio_ids,
                                           angio_symbols, cohort_dir, cohort_label) {
  angio_detected <- intersect(angio_ids, rownames(vsd_mat))
  if (length(angio_detected) == 0) {
    writeLines("Angiogenesis PSMA analysis skipped: no angiogenesis genes in VST matrix.",
               file.path(cohort_dir, paste0(cohort_label, "_angio_psma_summary.txt")))
    return(invisible(NULL))
  }
  if (!"PSMA_group" %in% colnames(cohort_meta) ||
      length(unique(na.omit(cohort_meta$PSMA_group))) < 2) {
    writeLines("Angiogenesis PSMA analysis skipped: PSMA_group absent or only 1 level.",
               file.path(cohort_dir, paste0(cohort_label, "_angio_psma_summary.txt")))
    return(invisible(NULL))
  }
  
  common_samples <- intersect(cohort_meta$sample_barcode, colnames(vsd_mat))
  z_raw <- row_zscore(vsd_mat[angio_detected, common_samples, drop = FALSE])
  angio_score <- colMeans(z_raw, na.rm = TRUE)
  cohort_meta$Angio_score <- angio_score[cohort_meta$sample_barcode]
  
  write.csv(cohort_meta %>%
              select(sample_barcode, patient_id, sample_type,
                     PSMA_group, FOLH1_vst, Angio_score),
            file.path(cohort_dir, paste0(cohort_label, "_angio_psma_scores.csv")),
            row.names = FALSE)
  
  make_signature_heatmap_clean(
    vsd_mat = vsd_mat,
    cohort_meta = cohort_meta,
    signature_ids = angio_ids,
    signature_symbols = angio_symbols,
    cohort_dir = cohort_dir,
    cohort_label = cohort_label,
    output_name = paste0(cohort_label, "_angio_psma_top25_heatmap"),
    signature_label = "Hallmark Angiogenesis",
    group_var = "PSMA_group",
    split_label = "PSMA Low"
  )
  
  cohort_meta$PSMA_group <- factor(cohort_meta$PSMA_group,
                                   levels = c("PSMA Low", "PSMA High"))
  wt <- safe_wilcox(Angio_score ~ PSMA_group, data = cohort_meta)
  
  p_box <- ggplot(cohort_meta,
                  aes(x = PSMA_group, y = Angio_score, fill = PSMA_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.55) +
    scale_fill_manual(values = c("PSMA Low" = "#4393C3", "PSMA High" = "#D6604D")) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = 7, color = "grey40"),
          plot.margin = margin(t = 20, r = 10, b = 10, l = 10)) +
    labs(title = paste0(cohort_label,
                        " - Hallmark Angiogenesis score (PSMA High vs Low)"),
         subtitle = paste0("HALLMARK_ANGIOGENESIS - ", length(angio_detected),
                           " genes measured (of ", length(angio_symbols), " in set)"),
         x = "", y = "Angiogenesis score (mean z-score)")
  p_box <- add_ttest_annotation(p_box, cohort_meta, "Angio_score",
                                "PSMA_group", "PSMA Low", "PSMA High")
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_angio_psma_boxplot.png")),
            p_box, width = 6.5, height = 5.5)
  
  p_scatter <- ggplot(cohort_meta,
                      aes(x = FOLH1_vst, y = Angio_score, color = PSMA_group)) +
    geom_point(size = 2.1, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7,
                linetype = "dashed", show.legend = FALSE) +
    scale_color_manual(values = c("PSMA Low" = "#4393C3", "PSMA High" = "#D6604D")) +
    theme_classic(base_size = 13) +
    labs(title = paste0(cohort_label, " - FOLH1 vs Angiogenesis score"),
         x = "FOLH1 VST expression",
         y = "Angiogenesis score (mean z-score)",
         color = "")
  save_plot(file.path(cohort_dir, paste0(cohort_label, "_FOLH1_vs_angio_scatter.png")),
            p_scatter, width = 5.8, height = 4.5)
  
  write.csv(data.frame(
    cohort = cohort_label,
    n_angio_genes = length(angio_detected),
    n_psma_low = sum(cohort_meta$PSMA_group == "PSMA Low", na.rm = TRUE),
    n_psma_high = sum(cohort_meta$PSMA_group == "PSMA High", na.rm = TRUE),
    median_low = median(cohort_meta$Angio_score[cohort_meta$PSMA_group == "PSMA Low"],
                        na.rm = TRUE),
    median_high = median(cohort_meta$Angio_score[cohort_meta$PSMA_group == "PSMA High"],
                         na.rm = TRUE),
    wilcox_p = if (is.null(wt)) NA_real_ else wt$p.value),
    file.path(cohort_dir, paste0(cohort_label, "_angio_psma_stats.csv")),
    row.names = FALSE)
  
  sym_map <- map_ensembl_to_annotations(angio_detected)
  write.csv(data.frame(
    ENSEMBL = angio_detected,
    SYMBOL = sym_map$SYMBOL[match(angio_detected, sym_map$ENSEMBL)]),
    file.path(cohort_dir, paste0(cohort_label, "_angio_psma_gene_list.csv")),
    row.names = FALSE)
}

run_angiogenesis_volcano <- function(res_shrunk_df, angio_ids, angio_symbols,
                                     cohort_dir, prefix) {
  if (is.null(res_shrunk_df) || nrow(res_shrunk_df) == 0 ||
      is.null(angio_ids) || length(angio_ids) == 0) {
    return(invisible(NULL))
  }
  
  vdf <- prepare_volcano_data_clean(res_shrunk_df,
                                    lfc_col = "log2FoldChange",
                                    padj_col = "padj",
                                    label_col = "SYMBOL",
                                    gene_id_col = "gene_id")
  vdf <- vdf %>%
    mutate(
      is_angio = gene_id %in% angio_ids,
      angio_class = dplyr::case_when(
        is_angio & direction == "up"   ~ "Angio - Up in Recurrent",
        is_angio & direction == "down" ~ "Angio - Up in Initial",
        is_angio                       ~ "Angio - NS",
        TRUE                           ~ "Other"
      ),
      angio_class = factor(angio_class,
                           levels = c("Angio - Up in Recurrent",
                                      "Angio - Up in Initial",
                                      "Angio - NS", "Other"))
    )
  
  n_det <- sum(vdf$is_angio, na.rm = TRUE)
  n_sig <- sum(vdf$is_angio & vdf$padj < clean_volcano_fdr_cutoff, na.rm = TRUE)
  if (n_det == 0) {
    writeLines("Angiogenesis volcano skipped: no genes detected.",
               file.path(cohort_dir, paste0(prefix, "_angio_volcano_summary.txt")))
    return(invisible(NULL))
  }
  
  write.csv(vdf %>% filter(is_angio),
            file.path(cohort_dir, paste0(prefix, "_angiogenesis_genes_DESeq2.csv")),
            row.names = FALSE)
  
  limits <- choose_volcano_limits_clean(vdf, "log2FoldChange")
  label_df <- vdf %>%
    filter(is_angio, direction != "ns") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 8)
  
  col_vals <- c("Angio - Up in Recurrent" = "#B2182B",
                "Angio - Up in Initial"   = "#2166AC",
                "Angio - NS"              = "#F4A582",
                "Other"                   = "grey82")
  
  p_a <- ggplot() +
    geom_point(data = vdf %>% filter(!is_angio),
               aes(log2FoldChange, neglog10),
               color = "grey78", alpha = 0.18, size = 0.75) +
    geom_point(data = vdf %>% filter(is_angio, direction == "ns"),
               aes(log2FoldChange, neglog10, color = angio_class),
               alpha = 0.82, size = 1.7) +
    geom_point(data = vdf %>% filter(is_angio, direction != "ns"),
               aes(log2FoldChange, neglog10, color = angio_class),
               alpha = 0.95, size = 2.1) +
    geom_vline(xintercept = c(-clean_volcano_lfc_cutoff, clean_volcano_lfc_cutoff),
               linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_hline(yintercept = -log10(clean_volcano_fdr_cutoff),
               linetype = "dashed", color = "grey45", linewidth = 0.35) +
    geom_text_repel(data = label_df,
                    aes(log2FoldChange, neglog10, label = label_text,
                        color = angio_class),
                    size = 3.0, min.segment.length = 0,
                    segment.size = 0.25, segment.alpha = 0.55,
                    box.padding = 0.35, point.padding = 0.25,
                    max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = col_vals, drop = FALSE) +
    coord_cartesian(xlim = limits$x, ylim = limits$y, clip = "off") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 15, face = "bold"),
          plot.subtitle = element_text(size = 8.8, color = "grey35"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 9.2),
          axis.title = element_text(size = 12.5),
          axis.text = element_text(size = 10.5),
          plot.margin = margin(t = 12, r = 18, b = 10, l = 10)) +
    labs(title = paste0(prefix, " - Angiogenesis gene signature"),
         subtitle = sprintf(
           "HALLMARK_ANGIOGENESIS: %d/%d genes detected, %d significant (FDR<0.05)",
           n_det, length(angio_ids), n_sig),
         x = "Shrunken log2 fold-change",
         y = "-log10 adjusted p-value")
  save_plot(file.path(cohort_dir, paste0(prefix, "_angiogenesis_volcano.png")),
            p_a, width = 7.4, height = 5.4)
  
  writeLines(c(paste("Prefix:", prefix),
               paste("Gene set size:", length(angio_ids)),
               paste("Detected:", n_det),
               paste("Sig FDR<0.05:", n_sig),
               paste("Up Recurrent:", sum(vdf$is_angio & vdf$direction == "up")),
               paste("Up Initial:", sum(vdf$is_angio & vdf$direction == "down"))),
             file.path(cohort_dir, paste0(prefix, "_angio_volcano_summary.txt")))
}

write_initial_vs_recurrent_fpkm_heatmap_clean <- function(counts_mat, counts_meta,
                                                          fpkm_mat, fpkm_meta,
                                                          gene_map_fpkm,
                                                          psm_result = NULL) {
  if (is.null(fpkm_mat) || is.null(fpkm_meta)) return(invisible(NULL))
  
  cohort_dir <- file.path(output_root, "Initial_vs_Recurrent")
  if (!is.null(psm_result)) {
    ini_bc <- psm_result$initial_barcodes
    rec_bc <- psm_result$recurrent_barcodes
  } else {
    ini_bc <- counts_meta %>% filter(sample_type_code == "01") %>% pull(sample_barcode)
    rec_bc <- counts_meta %>% filter(sample_type_code == "02") %>% pull(sample_barcode)
  }
  comp_samples <- c(ini_bc, rec_bc)
  if (length(comp_samples) < 2) return(invisible(NULL))
  
  counts_sub <- counts_mat[, comp_samples, drop = FALSE]
  counts_sub <- counts_sub[rowSums(counts_sub) > 0, , drop = FALSE]
  dds_norm <- DESeqDataSetFromMatrix(countData = round(counts_sub),
                                     colData = data.frame(row.names = comp_samples),
                                     design = ~ 1)
  dds_norm <- dds_norm[rowSums(counts(dds_norm)) > 10, ]
  dds_norm <- estimateSizeFactors(dds_norm)
  vsd_mat <- assay(vst(dds_norm, blind = TRUE))
  if (!folh1_id %in% rownames(vsd_mat)) return(invisible(NULL))
  
  comp_meta <- counts_meta %>%
    filter(sample_barcode %in% colnames(vsd_mat),
           sample_barcode %in% comp_samples) %>%
    arrange(match(sample_barcode, colnames(vsd_mat))) %>%
    mutate(Comparison_group = factor(sample_type, levels = c("Initial", "Recurrent")),
           FOLH1_vst = as.numeric(vsd_mat[folh1_id, sample_barcode]))
  
  common_fpkm <- intersect(comp_meta$sample_barcode, rownames(fpkm_mat))
  if (length(common_fpkm) < 2) return(invisible(NULL))
  
  meta_sub <- comp_meta %>%
    filter(sample_barcode %in% common_fpkm) %>%
    arrange(FOLH1_vst)
  fpkm_sub <- fpkm_mat[meta_sub$sample_barcode, , drop = FALSE]
  
  heat_mat <- t(log2(fpkm_sub + 1))
  heat_gene_ids <- rownames(heat_mat)
  heat_labels <- gene_map_fpkm$SYMBOL[match(heat_gene_ids, gene_map_fpkm$ENSEMBL)]
  rownames(heat_mat) <- make.unique(ifelse(is.na(heat_labels) | heat_labels == "",
                                           heat_gene_ids, heat_labels))
  heat_z <- row_zscore(heat_mat)
  ordered_s <- meta_sub$sample_barcode
  ordered_g <- cluster_order_from_z(heat_z[, ordered_s, drop = FALSE], "rows")
  hmap_df <- z_matrix_to_long(heat_z, sample_order = ordered_s,
                              gene_order = ordered_g)
  
  p_h <- ggplot(hmap_df, aes(sample_barcode, gene_label, fill = z)) +
    geom_tile() +
    z_heatmap_fill() +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          plot.subtitle = element_text(size = 8, color = "grey40")) +
    labs(title = "Selected genes (FPKM): Initial vs Recurrent",
         subtitle = "Samples ordered by FOLH1 expression low-to-high; genes clustered by z-score",
         x = paste0("Samples (n=", length(ordered_s), ")"),
         y = "Gene", fill = "z-score")
  save_plot(file.path(cohort_dir, "Initial_vs_Recurrent_fpkm_heatmap.png"),
            p_h, width = 8, height = 4.8)
}

# Wrap the original cohort function so the final pipeline also produces the new
# angiogenesis heatmap/score block for Initial and Recurrent cohorts.
if (!exists("analyse_cohort_original_for_clean_plots", inherits = FALSE)) {
  analyse_cohort_original_for_clean_plots <- analyse_cohort
}

analyse_cohort <- function(cohort_code, cohort_label,
                           counts_mat, counts_meta, fpkm_mat, fpkm_meta,
                           survival_df, hallmark_hypoxia_ids,
                           hallmark_hypoxia_symbols, gene_map_all, gene_map_fpkm,
                           angiogenesis_ids = NULL, angiogenesis_symbols = NULL) {
  analyse_cohort_original_for_clean_plots(
    cohort_code, cohort_label,
    counts_mat, counts_meta, fpkm_mat, fpkm_meta,
    survival_df, hallmark_hypoxia_ids, hallmark_hypoxia_symbols,
    gene_map_all, gene_map_fpkm
  )
  
  if (is.null(angiogenesis_ids) && exists("hallmark_angio_ids", inherits = TRUE)) {
    angiogenesis_ids <- get("hallmark_angio_ids", inherits = TRUE)
  }
  if (is.null(angiogenesis_symbols) && exists("hallmark_angio_symbols", inherits = TRUE)) {
    angiogenesis_symbols <- get("hallmark_angio_symbols", inherits = TRUE)
  }
  if (is.null(angiogenesis_ids) || length(angiogenesis_ids) == 0) {
    return(invisible(NULL))
  }
  
  cohort_dir <- file.path(output_root, cohort_label)
  cohort_samples <- counts_meta %>%
    filter(sample_type_code == cohort_code) %>%
    pull(sample_barcode)
  if (length(cohort_samples) < 2) return(invisible(NULL))
  
  counts_sub <- counts_mat[, cohort_samples, drop = FALSE]
  counts_sub <- counts_sub[rowSums(counts_sub) > 0, , drop = FALSE]
  dds_norm <- DESeqDataSetFromMatrix(countData = round(counts_sub),
                                     colData = data.frame(row.names = cohort_samples),
                                     design = ~ 1)
  dds_norm <- dds_norm[rowSums(counts(dds_norm)) > 10, ]
  dds_norm <- estimateSizeFactors(dds_norm)
  vsd_mat <- assay(vst(dds_norm, blind = TRUE))
  if (!folh1_id %in% rownames(vsd_mat)) return(invisible(NULL))
  
  cohort_meta <- counts_meta %>%
    filter(sample_barcode %in% colnames(vsd_mat)) %>%
    arrange(match(sample_barcode, colnames(vsd_mat))) %>%
    mutate(
      FOLH1_vst = as.numeric(vsd_mat[folh1_id, sample_barcode]),
      PSMA_group = make_psma_groups(FOLH1_vst)
    )
  rownames(cohort_meta) <- cohort_meta$sample_barcode
  
  run_angiogenesis_psma_analysis(vsd_mat, cohort_meta,
                                 angiogenesis_ids, angiogenesis_symbols,
                                 cohort_dir, cohort_label)
}

if (!exists("run_initial_vs_recurrent_analysis_original_for_clean_plots", inherits = FALSE)) {
  run_initial_vs_recurrent_analysis_original_for_clean_plots <- run_initial_vs_recurrent_analysis
}

run_initial_vs_recurrent_analysis <- function(
    counts_mat, counts_meta, fpkm_mat, fpkm_meta,
    survival_df, hallmark_hypoxia_ids, hallmark_hypoxia_symbols,
    angiogenesis_ids, angiogenesis_symbols,
    gene_map_all, gene_map_fpkm,
    psm_result = NULL) {
  run_initial_vs_recurrent_analysis_original_for_clean_plots(
    counts_mat = counts_mat,
    counts_meta = counts_meta,
    fpkm_mat = fpkm_mat,
    fpkm_meta = fpkm_meta,
    survival_df = survival_df,
    hallmark_hypoxia_ids = hallmark_hypoxia_ids,
    hallmark_hypoxia_symbols = hallmark_hypoxia_symbols,
    angiogenesis_ids = angiogenesis_ids,
    angiogenesis_symbols = angiogenesis_symbols,
    gene_map_all = gene_map_all,
    gene_map_fpkm = gene_map_fpkm,
    psm_result = psm_result
  )
  
  write_initial_vs_recurrent_fpkm_heatmap_clean(
    counts_mat = counts_mat,
    counts_meta = counts_meta,
    fpkm_mat = fpkm_mat,
    fpkm_meta = fpkm_meta,
    gene_map_fpkm = gene_map_fpkm,
    psm_result = psm_result
  )
}

# -----------------------------------------------------------------------------
# 15/16. Data loading (skipped if objects are already in memory)
# -----------------------------------------------------------------------------

if (!exists("clinical_df")) {
  message("Loading clinical data...")
  clinical_df <- read.csv(clinical_file, check.names = FALSE)
  survival_df <- prepare_survival_df(clinical_df)
} else {
  message("clinical_df already in environment - skipping.")
  if (!exists("survival_df"))
    survival_df <- prepare_survival_df(clinical_df)
}

if (!exists("counts_mat")) {
  message("Loading raw counts (this takes a while)...")
  counts_raw   <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  counts_dedup <- deduplicate_by_patient_and_type(counts_raw, rownames(counts_raw))
  counts_meta  <- counts_dedup$meta
  counts_mat   <- t(as.matrix(counts_dedup$data))
  rownames(counts_mat) <- clean_ensembl(rownames(counts_mat))
  storage.mode(counts_mat) <- "double"
  counts_mat[is.na(counts_mat)] <- 0
} else {
  message("counts_mat already in environment - skipping.")
}

if (!exists("fpkm_mat")) {
  message("Loading FPKM...")
  fpkm_raw   <- read.csv(fpkm_file, row.names = 1, check.names = FALSE)
  fpkm_dedup <- deduplicate_by_patient_and_type(fpkm_raw, rownames(fpkm_raw))
  fpkm_meta  <- fpkm_dedup$meta
  fpkm_mat   <- as.matrix(fpkm_dedup$data)
  storage.mode(fpkm_mat) <- "double"
  colnames(fpkm_mat) <- clean_ensembl(colnames(fpkm_mat))
} else {
  message("fpkm_mat already in environment - skipping.")
}

if (!exists("gene_map_all")) {
  message("Building gene annotations...")
  gene_map_all  <- map_ensembl_to_annotations(rownames(counts_mat))
  gene_map_fpkm <- map_ensembl_to_annotations(colnames(fpkm_mat))
} else {
  message("gene_map_all already in environment - skipping.")
}

if (!exists("hallmark_hypoxia_ids")) {
  message("Fetching MSigDB gene sets...")
  hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H")
  hallmark_hypoxia_symbols <- hallmark_sets %>%
    filter(gs_name == "HALLMARK_HYPOXIA") %>% pull(gene_symbol) %>% unique()
  hallmark_hypoxia_ids <- gene_map_all %>%
    filter(SYMBOL %in% hallmark_hypoxia_symbols) %>% pull(ENSEMBL) %>% unique()
  hallmark_angio_symbols <- hallmark_sets %>%
    filter(gs_name == "HALLMARK_ANGIOGENESIS") %>% pull(gene_symbol) %>% unique()
  hallmark_angio_ids <- gene_map_all %>%
    filter(SYMBOL %in% hallmark_angio_symbols) %>% pull(ENSEMBL) %>% unique()
} else {
  message("Hallmark gene sets already in environment - skipping.")
}

if (!exists("clinical_psm_df")) {
  message("Preparing PSM covariates...")
  clinical_psm_df <- prepare_clinical_for_psm(clinical_df)
} else {
  message("clinical_psm_df already in environment - skipping.")
}

# -----------------------------------------------------------------------------
# 17. Run analyses
# -----------------------------------------------------------------------------

message("=== Analysing Initial cohort ===")
analyse_cohort("01", "Initial", counts_mat, counts_meta, fpkm_mat, fpkm_meta,
               survival_df, hallmark_hypoxia_ids, hallmark_hypoxia_symbols,
               gene_map_all, gene_map_fpkm)

message("=== Analysing Recurrent cohort ===")
analyse_cohort("02", "Recurrent", counts_mat, counts_meta, fpkm_mat, fpkm_meta,
               survival_df, hallmark_hypoxia_ids, hallmark_hypoxia_symbols,
               gene_map_all, gene_map_fpkm)

# ── TIER A : Paired analysis on the 9 patients with both samples ─────────────
message("=== TIER A: Paired DESeq2 analysis (~ patient_id + sample_type) ===")
run_paired_deseq_analysis(
  counts_mat           = counts_mat,
  counts_meta          = counts_meta,
  hallmark_hypoxia_ids = hallmark_hypoxia_ids,
  hallmark_hypoxia_symbols = hallmark_hypoxia_symbols,
  angiogenesis_ids     = hallmark_angio_ids,
  angiogenesis_symbols = hallmark_angio_symbols,
  gene_map_all         = gene_map_all)

# ── TIER B : PSM 1:2 matched cohort (13 recurrent vs 26 matched initial) ─────
message("=== TIER B: Propensity Score Matching (1:2 — 2 initial per recurrent) ===")
psm_result <- run_propensity_matching(
  counts_meta, clinical_psm_df,
  psm_dir = file.path(output_root, "PSM_1to2"),
  ratio   = psm_ratio)   # psm_ratio = 2L

message("=== TIER B: Initial vs Recurrent matched cohort analysis ===")
run_initial_vs_recurrent_analysis(
  counts_mat               = counts_mat,
  counts_meta              = counts_meta,
  fpkm_mat                 = fpkm_mat,
  fpkm_meta                = fpkm_meta,
  survival_df              = survival_df,
  hallmark_hypoxia_ids     = hallmark_hypoxia_ids,
  hallmark_hypoxia_symbols = hallmark_hypoxia_symbols,
  angiogenesis_ids         = hallmark_angio_ids,
  angiogenesis_symbols     = hallmark_angio_symbols,
  gene_map_all             = gene_map_all,
  gene_map_fpkm            = gene_map_fpkm,
  psm_result               = psm_result)

message("=== DONE. Outputs in: ", output_root, " ===")
message("  Tier A (paired):  ", file.path(output_root, "Paired_Analysis"))
message("  Tier B (matched): ", file.path(output_root, "Initial_vs_Recurrent"))
message("  PSM balance:      ", file.path(output_root, "PSM_1to2"))
