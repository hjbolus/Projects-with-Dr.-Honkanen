# load GSEA utils
if (!exists("gsea_utils_initialized", envir=globalenv())) {source("~/path/to/GSEA utils.r")}

# initialize output directory if needed
output_dir <- "output/directory/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# load nanopore data (ensembl gene IDs x gene counts for each sample)
if (!exists("nanopore_results_table", envir=globalenv())) {
  nanopore_df <- read_excel(".../Nanopore DE analysis.xlsx", 
                         sheet = "...", 
                         range = "...",
                         col_types = c(...)
                         )
  
  # strip ensembl ID version numbers
  nanopore_df$gene_id <- sub("\\..*$", "", nanopore_df$gene_id)
  
  # select columns containing count data
  nanopore_counts_matrix <- nanopore_df %>%
    dplyr::select(ends_with("_count")) %>%
    as.matrix()
  
  # use edgeR filterbyexpr() to remove low expression genes
  rownames(nanopore_counts_matrix) <- nanopore_df$gene_id
  nanopore_groups <- factor(c(rep("CONTROL", 4), rep("VARIANT", 4)))
  nanopore_dge <- DGEList(
    counts = nanopore_counts_matrix,
    genes = nanopore_df[, c("gene_id", "gene_name")]  # keep annotation
  )
  
  keep <- filterByExpr(nanopore_dge, group = nanopore_groups)
  nanopore_dge <- nanopore_dge[keep, , keep.lib.sizes = FALSE]
  
  # recalculate normalized FC and p-value after filtering
  nanopore_dge <- calcNormFactors(nanopore_dge)
  nanopore_dge$samples$group <- nanopore_groups
  nanopore_design <- model.matrix(~ nanopore_groups)
  
  nanopore_dge <- estimateDisp(nanopore_dge, nanopore_design)
  
  nanopore_fit <- glmQLFit(nanopore_dge, nanopore_design)
  nanopore_qlf <- glmQLFTest(nanopore_fit, coef = 2)  # coef=2 corresponds to group VARIANT
  
  nanopore_results <- topTags(nanopore_qlf, n = Inf)  # all genes
  nanopore_results_table <- nanopore_results$table    # data.frame with logFC, logCPM, F, PValue, FDR
}

# calculate ranks
c(nanopore_results_table, nanopore_ranks_signed_p) %<-% calc_ranks(
  df = nanopore_results_table,
  id_col = 'gene_id',
  logfc_col = 'logFC', 
  p_col = 'PValue',
  rank_type = 'signed_p'
)

# plot ranks against log2FC and -log10(p-value). Available metrics are "signed_p", "pi_stat", and "signed_p".
save_gsea_plot(plot_metric(nanopore_results_table, "signed_p"), paste0(output_dir, "nanopore_ranks_signed_p.png"))

# run fGSEA
nanopore_signed_p_h <- run_and_plot_fgsea(h_gs_gsymbol, nanopore_ranks_signed_p)

# save figures
save_sig_enrichment_plots(results = nanopore_signed_p_h,
                          genesets = h_gs_gsymbol,
                          ranks = nanopore_ranks_signed_p,
                          path = output_dir,
                          name = "nanopore_signed_p_h")

# save complete results
save(nanopore_signed_p_h, file=paste0(output_dir, "nanopore_signed_p_h", ".RData"))
