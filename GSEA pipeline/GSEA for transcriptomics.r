needed <- c("msigdbr", "fgsea", "dplyr", "edgeR", "ggplot2")
to_install <- needed[!needed %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

library(msigdbr)
library(fgsea)
library(dplyr)
library(readxl)
library(edgeR)
library(ggplot2)
source("... GSEA utils.r")

# declare output path
output_path = "..."

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

# calculate and assign ranks
nanopore_results_table$pi_score <- nanopore_results_table$logFC * -log(nanopore_results_table$PValue, 10)
nanopore_ranks <- nanopore_results_table$pi_score
names(nanopore_ranks) <- nanopore_results_table$gene_id
nanopore_ranks <- nanopore_ranks + rank(names(nanopore_ranks)) * 1e-12
sort(nanopore_ranks, decreasing=TRUE)

# plot rank metric against log2FC and -log10(p-value)
save_gsea_plot(plot_metric(nanopore_results_table, "pi_score"), paste(output_path, "nanopore pi score.png", sep=''))

# run fgsea of selected gene sets
fgsea_nano_tft <- run_and_plot_fgsea(tft_gs, nanopore_ranks)
fgsea_nano_h <- run_and_plot_fgsea(h_gs, nanopore_ranks)
fgsea_nano_mirna <- run_and_plot_fgsea(mirna_gs, nanopore_ranks)
fgsea_nano_rea <- run_and_plot_fgsea(reactome_gs, nanopore_ranks)
fgsea_nano_k <- run_and_plot_fgsea(kegg_gs, nanopore_ranks)
fgsea_nano_gobp <- run_and_plot_fgsea(gobp_gs, nanopore_ranks)
fgsea_nano_gomf <- run_and_plot_fgsea(gomf_gs, nanopore_ranks)
fgsea_nano_gocc <- run_and_plot_fgsea(gocc_gs, nanopore_ranks)

# save results
save_sig_enrichment_plots(results = fgsea_nano_tft, 
                          genesets = tft_gs, 
                          ranks = nanopore_ranks, 
                          path = output_path, 
                          name = "fgsea_nano_tft")
                          
save_sig_enrichment_plots(results = fgsea_nano_h,
                          genesets = h_gs,
                          ranks = nanopore_ranks,
                          path = output_path,
                          name = "fgsea_nano_h")

save_sig_enrichment_plots(results = fgsea_nano_mirna,
                          genesets = mirna_gs,
                          ranks = nanopore_ranks,
                          path = output_path,
                          name = "fgsea_nano_mirna")

save_sig_enrichment_plots(results = fgsea_nano_rea,
                          genesets = reactome_gs,
                          ranks = nanopore_ranks,
                          path = output_path,
                          name = "fgsea_nano_rea")

save_sig_enrichment_plots(results = fgsea_nano_k,
                          genesets = kegg_gs,
                          ranks = nanopore_ranks,
                          path = output_path,
                          name = "fgsea_nano_k")


save_sig_enrichment_plots(results = fgsea_nano_gobp,
                          genesets = gobp_gs,
                          ranks = nanopore_ranks,
                          path = output_path,
                          name = "fgsea_nano_gobp")

save_sig_enrichment_plots(results = fgsea_nano_gomf,
                          genesets = gomf_gs,
                          ranks = nanopore_ranks,
                          path = output_path,
                          name = "fgsea_nano_gomf")

save_sig_enrichment_plots(results = fgsea_nano_gocc,
                          genesets = gocc_gs,
                          ranks = nanopore_ranks,
                          path = output_path,
                          name = "fgsea_nano_gocc")

