# load GSEA utils
if (!exists("gsea_utils_initialized", envir=globalenv())) {source("~/path/to/GSEA utils.r")}

# initialize output directory if needed
output_dir <- "output/directory/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# load prot data
if (!exists("prot_results_table", envir=globalenv())) {
  prot_results_table <- "load proteomics results table with gene ids as Gene Symbols or Ensembl IDs, p-values, and log2 FC"
  
  colnames(heke198k_prot_results_table) <- c('gene_id', 'PValue', 'logFC')
}

# calculate ranks
c(prot_results_table, prot_ranks_pi_stat) %<-% calc_ranks(
  df = prot_results_table,
  id_col = 'gene_id',
  logfc_col = 'logFC', 
  p_col = 'PValue',
  rank_type = 'pi_stat'
)

# plot ranks against log2FC and -log10(p-value). Available metrics are "signed_p", "pi_stat", and "msd".
save_gsea_plot(plot_metric(prot_results_table, "pi_stat"), paste0(output_dir, "prot_ranks_pi_stat.png"))

# run fGSEA
prot_signed_p_h <- run_and_plot_fgsea(h_gs_gsymbol, prot_ranks_pi_stat)

# save figures
save_sig_enrichment_plots(results = prot_pi_stat_h,
                          genesets = h_gs_gsymbol,
                          ranks = prot_ranks_pi_stat,
                          path = output_dir,
                          name = "prot_pi_stat_h")

# save complete results
save(prot_pi_stat_h, file=paste0(output_dir, "prot_pi_stat_h", ".RData"))
