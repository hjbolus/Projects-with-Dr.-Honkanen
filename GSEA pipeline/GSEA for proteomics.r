needed <- c("msigdbr", "fgsea", "dplyr", "edgeR", "ggplot2")
to_install <- needed[!needed %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

library(msigdbr)
library(fgsea)
library(dplyr)
library(readxl)
library(ggplot2)
source("... GSEA utils")


# declare output path
output_path = "..."

# load proteomics data - only gene name, p-value, and logFC needed
if (!exists("prot_results_table", envir=globalenv())) {
  prot_results_table <- read_excel(".../proteomics.xlsx", 
                                   sheet = "...", 
                                   range = "...", 
                                  col_types = c(...)
                                  )
  
  colnames(prot_results_table) <- c('gene id', 'PValue', 'logFC')
}

# calculate ranks
prot_results_table$pi_score <- prot_results_table$logFC * -log(prot_results_table$PValue, 10)
prot_ranks_pi <- prot_results_table$pi_score
names(prot_ranks_pi) <- prot_results_table$`gene id`
prot_ranks_pi <- prot_ranks_pi[!is.na(prot_ranks_pi)]
prot_ranks_pi <- prot_ranks_pi + rank(names(prot_ranks_pi)) * 1e-12 # break ties alphabetically
prot_ranks_pi <- tapply(prot_ranks_pi, names(prot_ranks_pi), function(x) x[which.max(abs(x))])
sort(prot_ranks_pi, decreasing=TRUE)

# plot rank metrics against log2FC and -log10(p-value)
save_gsea_plot(plot_metric(prot_results_table, "pi_score"), paste(output_path, "prot pi score.png", sep='')

# run fgsea of selected genesets
fgsea_prot_tft_pi <- run_and_plot_fgsea(tft_gs_gsymbol, prot_ranks_pi)
fgsea_prot_h_pi <- run_and_plot_fgsea(h_gs_gsymbol, prot_ranks_pi)
fgsea_prot_mirna_pi <- run_and_plot_fgsea(mirna_gs_gsymbol, prot_ranks_pi)
fgsea_prot_rea_pi <- run_and_plot_fgsea(reactome_gs_gsymbol, prot_ranks_pi)
fgsea_prot_k_pi <- run_and_plot_fgsea(kegg_gs_gsymbol, prot_ranks_pi)
fgsea_prot_gobp_pi <- run_and_plot_fgsea(gobp_gs_gsymbol, prot_ranks_pi)
fgsea_prot_gomf_pi <- run_and_plot_fgsea(gomf_gs_gsymbol, prot_ranks_pi)
fgsea_prot_gocc_pi <- run_and_plot_fgsea(gocc_gs_gsymbol, prot_ranks_pi)

# save results
save_sig_enrichment_plots(results = fgsea_prot_tft_pi,
                          genesets = tft_gs_gsymbol,
                          ranks = prot_ranks_pi,
                          path = output_path,
                          name = "fgsea_prot_tft")

save_sig_enrichment_plots(results = fgsea_prot_h_pi,
                          genesets = h_gs_gsymbol,
                          ranks = prot_ranks_pi,
                          path = output_path,
                          name = "fgsea_prot_h")

save_sig_enrichment_plots(results = fgsea_prot_mirna_pi,
                          genesets = mirna_gs_gsymbol,
                          ranks = prot_ranks_pi,
                          path = output_path,
                          name = "fgsea_prot_mirna")

save_sig_enrichment_plots(results = fgsea_prot_rea_pi,
                          genesets = reactome_gs_gsymbol,
                          ranks = prot_ranks_pi,
                          path = output_path,
                          name = "fgsea_prot_rea")

save_sig_enrichment_plots(results = fgsea_prot_k_pi,
                          genesets = kegg_gs_gsymbol,
                          ranks = prot_ranks_pi,
                          path = output_path,
                          name = "fgsea_prot_k")

save_sig_enrichment_plots(results = fgsea_prot_gobp_pi,
                          genesets = gobp_gs_gsymbol,
                          ranks = prot_ranks_pi,
                          path = output_path,
                          name = "fgsea_prot_gobp")

save_sig_enrichment_plots(results = fgsea_prot_gomf_pi,
                          genesets = gomf_gs_gsymbol,
                          ranks = prot_ranks_pi,
                          path = output_path,
                          name = "fgsea_prot_gomf")

save_sig_enrichment_plots(results = fgsea_prot_gocc_pi,
                          genesets = gocc_gs_gsymbol,
                          ranks = prot_ranks_pi,
                          path = output_path,
                          name = "fgsea_prot_gocc")
