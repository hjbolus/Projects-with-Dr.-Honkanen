# load GSEA utils
if (!exists("gsea_utils_initialized", envir=globalenv())) {source("~/path/to/GSEA utils.r")}

# initialize output directory if needed
output_dir <- "output/directory/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# load phos data
if (!exists("phos_results_table", envir=globalenv())) {
  phos_results_table <- "load phosphoproteomics results table with gene ids as Gene Symbols or Ensembl IDs, p-values, and log2 FC"
  
  colnames(heke198k_phos_results_table) <- c('gene_id', 'PValue', 'logFC')
}

# calculate ranks
c(phos_results_table, phos_ranks_pi_stat) %<-% calc_ranks(
  df = phos_results_table,
  id_col = 'gene_id',
  logfc_col = 'logFC', 
  p_col = 'PValue',
  rank_type = 'pi_stat',
  abs = TRUE
)

# plot ranks against log2FC and -log10(p-value). Available metrics are "signed_p", "pi_stat", and "msd".
save_gsea_plot(plot_metric(phos_results_table, "pi_stat"), paste0(output_dir, "phos_ranks_pi_stat.png"))

# run fGSEA. scoreType = "pos" is important for phosphoproteomics, where enrichment/depletion does not directly map onto activity
phos_signed_p_h <- run_and_plot_fgsea(h_gs_gsymbol, phos_ranks_pi_stat, scoreType="pos")

# save figures
save_sig_enrichment_plots(results = phos_pi_stat_h,
                          genesets = h_gs_gsymbol,
                          ranks = phos_ranks_pi_stat,
                          path = output_dir,
                          name = "phos_pi_stat_h")

# save complete results
save(phos_signed_p_h, file=paste0(output_dir, "phos_pi_stat_h", ".RData"))
