# helper functions and statements for HEK E198K GSEA pipeline

needed <- c("msigdbr", "fgsea", "dplyr", "edgeR", "ggplot2", "stringr", "poolr")
to_install <- needed[!needed %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
# fgsea/msigdbr are on CRAN; if you prefer Bioc, you can install from Bioconductor too.

library(msigdbr)
library(fgsea)
library(dplyr)
library(readxl)
library(edgeR)
library(ggplot2)
library(stringr)
library(poolr)
library(writexl)
library(biomaRt)
library(WebGestaltR)
library(ggrepel)

# When plotting volcano, include independent pathways among top 10 positive and top 10 negative NES and add an asterisk in front of independent pathways: 362

biomart_table <- read_excel("~/Desktop/Research/Dr. Honkanen/Reference tables/biomart table 012625.xlsx", 
                                        col_types = c("text", "skip", "skip", 
                                                      "skip", "skip", "skip", "text", "skip", 
                                                      "skip", "skip", "skip"))

wrap.it <- function(x, len)
{ 
  sapply(x, function(y) paste(strwrap(y, len), 
                              collapse = "\n"), 
         USE.NAMES = FALSE)
}

is_listcolumn_all_ensembl <- function(gene_listcolumn) {
  # Regex for Ensembl gene IDs (with optional version number)
  ens_rx <- "^ENSG[0-9]{11}(\\.[0-9]+)?$"
  
  all(
    vapply(
      unlist(gene_listcolumn),             # flatten into character vector
      function(x) grepl(ens_rx, x),
      logical(1)
    )
  )
}

map_ensembl_listcolumn <- function(id_list, ref_table) {
  # ref_table must have columns: "Gene stable ID", "Gene name"
  
  # Build list: Ensembl ID -> vector of gene names (can be length > 1)
  id_to_names <- split(ref_table$`Gene name`, ref_table$`Gene stable ID`)
  
  lapply(id_list, function(vec) {
    # Strip version suffix, e.g. ENSG00000123456.7 → ENSG00000123456
    vec_clean <- sub("\\..*$", "", vec)
    
    # For each Ensembl ID, get all gene names (or the ID itself if no match)
    mapped_list <- lapply(seq_along(vec_clean), function(i) {
      id_clean <- vec_clean[i]
      original <- vec[i]
      
      if (!is.null(id_to_names[[id_clean]])) {
        unique(id_to_names[[id_clean]])  # all names for that Ensembl ID
      } else {
        original                         # keep Ensembl ID if no mapping
      }
    })
    
    # Flatten to a simple character vector
    unlist(mapped_list, use.names = FALSE)
  })
}

map_genename_listcolumn <- function(name_list, ref_table) {
  # Build lookup: gene name → ensembl ID
  lookup <- setNames(ref_table$`Gene stable ID`, ref_table$`Gene name`)
  
  # Apply mapping to each vector in list-column
  lapply(name_list, function(vec) {
    mapped <- lookup[vec]  # returns NA if no match
    
    # keep original name if no match
    mapped[is.na(mapped)] <- vec[is.na(mapped)]
    
    unname(mapped)
  })
}

convert_listcols_to_char <- function(df) {
  # if it's a data.table or tibble, that's fine; we just treat it like a data.frame
  listcols <- vapply(df, is.list, logical(1))
  listnames <- names(df)[listcols]
  
  if (!any(listcols)) return(df)
  
  # function to convert a single list-column to character
  conv_fun <- function(col) {
    vapply(col, function(x) {
      if (is.null(x)) {
        NA_character_
      } else if (length(x) == 0) {
        ""
      } else {
        toString(x)  # collapse vector: c("A","B") -> "A, B"
      }
    }, character(1))
  }
  
  # update each list-column in-place, column by column
  for (nm in listnames) {
    df[[nm]] <- conv_fun(df[[nm]])
  }
  
  df
}

convert_to_genesets <- function(collection, id_type='ensembl') {
  #id_type should be "ensembl" or "gene symbol". default ensembl
  prefixes <- c("GOBP_", "GOMF_", "GOCC_", "KEGG_", "REACTOME_", "HALLMARK_")
  pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
  
  id_type <- ifelse(id_type == 'ensembl', 'ensembl_gene',
                    ifelse(id_type == 'gene symbol', 'gene_symbol', FALSE))
  
  collection <- split(collection[[id_type]], collection$gs_name)

  new_names <- names(collection) |>
    sub(pattern, "", x = _) |>    # remove prefix
    gsub("_", " ", x = _)         # convert underscores to spaces
  
  # assign new names and return list
  names(collection) <- new_names
  return(collection)
}

combine_and_filter_genesets <- function(gs1, gs2) {
  return(
    rbind(gs1, gs2) %>%
    filter(!is.na(ensembl_gene) & ensembl_gene != "" & !grepl("UNKNOWN", gs_name)) %>%     # keep only rows with Ensembl IDs present, with a known TF
    mutate(ensembl_gene = sub("\\.\\d+$", "", ensembl_gene)) %>%                # strip version suffixes like ENSG00000123456.12
    mutate(gs_name = str_remove(gs_name, "_TARGET_GENES")) %>%                  #replace _TARGET_GENES with ''
    mutate(gs_name = as.character(gs_name)) %>%                                 # strip first suffix
    mutate(gs_name = sapply(strsplit(gs_name, "_"), function(x) {
        if (length(x) > 1) paste(head(x, -1), collapse = "_") else x
      })) %>%
    mutate(gs_name = str_remove(gs_name, "^[ACGTURYKMSWBDHVN]+_")) %>%          # remove motif names anchored to start, followed by an underscore
    mutate(gs_name = str_remove(gs_name, "_..$"))
  )
}

crit_val <- qnorm(0.975)

calc_msd <- function(logFC, p) {
  # based on https://www.bmj.com/content/343/bmj.d2090.extract
  # MSD is the minimum value of the 95% confidence interval
  # calculate z-score
  z <- qnorm(1 - p/2)
  
  # calculate standard error
  se <- abs(logFC / z)
  
  # calculate CI
  ci_left <- logFC + crit_val * se
  ci_right <- logFC - crit_val * se
  
  return(ifelse(abs(ci_left) < abs(ci_right), ci_left, ci_right))
}

plot_metric <- function(table, metric, x_max = 3, arrow_offset = 0.02, 
                        segment_length = 0.05, arrow_length = 0.1) {
  # x_max: absolute max for logFC to display
  # arrow_width: horizontal length of arrows (in data units)
  
  # Clamp x values and flag extreme points
  table2 <- table |>
    dplyr::mutate(
      xval       = logFC,
      x_clamped  = pmax(pmin(xval,  x_max), -x_max),
      is_clamped = abs(xval) > x_max
    )
  
  p <- ggplot(table2, aes(x = x_clamped, y = -log10(PValue))) +
    geom_point(
      aes(color = .data[[metric]], size = abs(.data[[metric]])),
      alpha = 0.8
    ) +
    
    # --- Arrows for points beyond +x_max ---
    geom_segment(
      data = subset(table2, xval > x_max),
      aes(
        x    =  x_max - arrow_offset - segment_length,
        xend =  x_max - arrow_offset,
        y    = -log10(PValue),
        yend = -log10(PValue)
      ),
      arrow = arrow(length = grid::unit(arrow_length, "cm"))
    ) +
    
    # Arrows for points beyond -x_max
    geom_segment(
      data = subset(table2, xval < -x_max),
      aes(
        x    = -x_max + arrow_offset + segment_length,
        xend = -x_max + arrow_offset,
        y    = -log10(PValue),
        yend = -log10(PValue)
      ),
      arrow = arrow(length = grid::unit(arrow_length, "cm"))
    ) +
    
    scale_x_continuous(
      limits = c(-x_max, x_max),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_color_gradient2(
      low = "blue", mid = "grey", high = "red", midpoint = 0,
      transform = "pseudo_log"
    ) +
    scale_size_continuous(range = c(1, 3)) +
    theme_minimal() +
    labs(
      x = "log2 Fold Change",
      y = "-log10(P-value)",
      color = metric,
      size = paste("|", metric, "|", sep = ""),
      title = "log2 FC vs -log10 pvalue"
    )
  
  return(p)
}


save_gsea_plot <- function(plot, name, height=1000, width=3000, res=300) {
  png(name, width, height, res = 150)
  print(plot)
  dev.off()
}

default_fgsea <- function(gs, ranks, scoreType = "std") {
  fgsea_res <- fgsea(
    pathways = gs,
    stats    = ranks,
    minSize  = 15,   # minimum genes in a set
    maxSize  = 500,  # maximum genes in a set
    scoreType = scoreType
    )
  fgsea_res <- fgsea_res[order(fgsea_res$padj), ]
  return(fgsea_res)
}

trim_gs <- function(gs, ranks, pct=0.2) {
  gs <- gs[sapply(gs, function(gs) {
    n_gs = length(gs)
    n_ranks = length(ranks)
    length(intersect(gs, names(ranks))) > max(10, ceiling(pct*n_gs))
  })]
  return(gs)
}

run_and_plot_fgsea <- function(genesets, ranks, scoreType = "std", pct=0.2) {
  genesets <- trim_gs(gs=genesets, ranks=ranks, pct=pct)
  fgsea_res <- default_fgsea(genesets, ranks, scoreType)
  
  if (is_listcolumn_all_ensembl(fgsea_res$leadingEdge)) {
    fgsea_res$leadingEdge.Ensembl <- fgsea_res$leadingEdge
    fgsea_res$leadingEdge.GeneName <- map_ensembl_listcolumn(fgsea_res$leadingEdge, biomart_table)
    fgsea_res$leadingEdge <- NULL
  } else {
    fgsea_res$leadingEdge.Ensembl <- map_genename_listcolumn(fgsea_res$leadingEdge, biomart_table)
    fgsea_res$leadingEdge.GeneName <- fgsea_res$leadingEdge
    fgsea_res$leadingEdge <- NULL
  }
  independent_pathways <- collapsePathways(fgseaRes = fgsea_res[fgsea_res$pval < 0.05],
                          stats = ranks,
                          pathways = genesets,
                          pval.threshold = 0.05
                          )
  ind <- fgsea_res[fgsea_res$pathway %in% independent_pathways$mainPathways]
  fgsea_res$ind <- fgsea_res$pathway %in% independent_pathways$mainPathways
  sig <- fgsea_res[fgsea_res$pval < 0.05]
  top10 <- fgsea_res[order(fgsea_res$padj, decreasing=FALSE), ][1:10, ]

  plot <- plotGseaTable(
    pathways = genesets[top10$pathway],
    stats    = ranks,
    fgseaRes = top10,
    gseaParam = 1
  )
  
  return(list(all=fgsea_res, sig=sig, ind=ind, top10=top10, plot=plot))
}

stouffer_method <- function(p) {
  poolr::stouffer(p)$p
}

wsc <- function(results, gs, pval_cutoff=0.05, topN=10) {
  
  pval_cols <- grep("^pval", names(results), value = TRUE)
  
  results$p.comb <- apply(results[, ..pval_cols], 1, stouffer_method)
  results$padj.comb <- results$padj.comb <- p.adjust(results$p.comb, method = "BH")
  
  fgsea_sig <- subset(results, padj.comb < 0.25) # this is just to determine which labels to display
  fgsea_sig <- fgsea_sig[order(fgsea_sig$padj.comb), ]
  set_names <- fgsea_sig$pathway
  idsInSet <- gs[set_names]
  pv <- fgsea_sig$padj.comb  # already in the same order as set_names
  pv[pv == 0] <- .Machine$double.xmin
  costs <- 1 / -log10(pv)
  names(costs) <- names(idsInSet)
  wsc_res <- weightedSetCover(
    idsInSet = idsInSet,
    costs    = costs,
    topN     = topN
  ) 
}

save_sig_enrichment_plots <- function(results, genesets, ranks, path, name) {
  newpath <- file.path(paste(path, name, sep=""))
  dir.create(newpath)

  df <- results$all
  
  #individual enrichment plots
  sig_pathways <- df[df$padj < 0.05]$pathway
  filtered_genesets <- genesets[names(genesets) %in% sig_pathways]
  for (gs in names(filtered_genesets)) {
    png(filename=paste(newpath, "/", gs, ".png", sep=""), height=750, width=1250, res=300)
    print(plotEnrichment(filtered_genesets[[gs]], ranks))
    dev.off()
  }
  
  # save excel of sig paths
  results$sig <- convert_listcols_to_char(results$sig)
  write_xlsx(results$sig, paste(newpath, "/", name, ".xlsx", sep=""))
  save(results, file=paste(newpath, "/", name, ".RData", sep=""))
  
  # gsea plot
  if (length(results$top10$pathway) > 0) {
  maxlen <- max(nchar(results$top10$pathway))
  width <- 3000
  if (!is.na(maxlen)){
    if (maxlen > 10) {
    width <- width + (maxlen-10) * 75
    }
  }
  
  png(filename=paste(newpath, "/gsea plot ", name, ".png", sep=""), res=300, height=1200, width=width)
  print(results$plot)
  dev.off()
  }
  
  df$pathway <- ifelse(
    isTRUE(df$ind) | (!is.na(df$ind) & df$ind),
    paste0("*", df$pathway),
    df$pathway
  )
  df$pathway <- wrap.it(df$pathway, 25)
  
  # volcano plot
  df$`-log10 p-value` <- -1 * log10(df$pval)
  df$`signed log2 NES` <- sign(df$NES)*log2(abs(df$NES))
  df$`Leading edge size` <- df$size
  df$`BH-adjusted p-value` <- df$padj
  df <- df %>% mutate(sig = padj <= 0.05)

  # initialize labels column
  df$label <- NA_character_
  
  # eligibility
  eligible <- !is.na(df$pval) & df$pval < 0.05 & !is.na(df$NES)
  
  pick_labels <- function(idx, n = 10) {
    if (!length(idx)) return(integer(0))
    
    # independent first
    idx1 <- idx[df$ind[idx]]
    idx1 <- idx1[order(df$pval[idx1], na.last = TRUE)]
    chosen <- head(idx1, n)
    
    # fill remainder by lowest p-values
    if (length(chosen) < n) {
      idx2 <- setdiff(idx, chosen)
      idx2 <- idx2[order(df$pval[idx2], na.last = TRUE)]
      chosen <- c(chosen, head(idx2, n - length(chosen)))
    }
    
    chosen
  }
  
  # positive NES
  pos_idx <- which(eligible & df$NES > 0)
  pos_labels <- pick_labels(pos_idx, 8)
  df$label[pos_labels] <- df$pathway[pos_labels]
  
  # negative NES
  neg_idx <- which(eligible & df$NES < 0)
  neg_labels <- pick_labels(neg_idx, 8)
  df$label[neg_labels] <- df$pathway[neg_labels]
  
  p <- ggplot() +
  # ----------------------
  # NON-SIGNIFICANT POINTS (padj >= 0.15)
  # ----------------------
  
  geom_point(
    data = subset(df, !sig),
    aes(
      x    = NES,#`signed log2 NES`,
      y    = `-log10 p-value`,
      color = `BH-adjusted p-value`,
      size  = `Leading edge size`
    ),
    alpha = 0.9
  ) +
    paletteer::scale_colour_paletteer_c(
      "grDevices::Purples 2",
      limits = c(0.05, 1),
      oob    = scales::squish,
      name   = "BH-adjusted p-value (> 0.05)"
    ) +
    
    ggnewscale::new_scale_color() +
    
    # ----------------------
  # SIGNIFICANT POINTS (padj < 0.05)
  # ----------------------
  
  geom_point(
    data = subset(df, sig),
    aes(fill = `BH-adjusted p-value`,, 
        size = `Leading edge size`,
        x    = NES,#`signed log2 NES`,
        y    = `-log10 p-value`),
    shape = 21,
    color = "#6D0026",
    stroke = 0.3,
    alpha = 0.9,
  ) +
    
    paletteer::scale_fill_paletteer_c(
      "grDevices::Reds",
      limits = c(0, 0.05),
      oob = scales::squish,
      name = "BH-adjusted p-value (≤ 0.05)"
    ) +
    
    # ----------------------
  # LABELS + HLINE + THEME
  # ----------------------
  geom_hline(
    yintercept = 1.30103,
    linetype   = "dashed",
    color      = "grey",
    linewidth  = 0.3,
  ) +
    theme_light()
  
  png(filename=paste(newpath, "/volcano plot ", name, ".png", sep=""), width=2000, height=1600, res=300)
  print(p)
  dev.off()
  
  p <- p +
    scale_x_continuous(expand = expansion(mult = 0.10)) +
    scale_y_continuous(expand = expansion(mult = 0.10)) +
    ggrepel::geom_label_repel(
      data = subset(df, !is.na(label)),
      aes(x = NES, y = `-log10 p-value`, label = label),
      
      # don't silently drop labels
      max.overlaps = Inf,
      
      # give the solver space + time
      max.time = 120,
      max.iter = 1e6,
      
      # increase separation
      box.padding   = 0.8,
      point.padding = 0.2,
      
      # push labels apart; reduce pull-back crowding
      force      = 1,
      force_pull = 1,
      
      # styling
      size = 1.2,
      alpha = 0.5,
      color = "black",
      min.segment.length = 0.1
    )
  
  png(filename=paste(newpath, "/volcano plot labeled ", name, ".png", sep=""), width=2000, height=1600, res=300)
  print(p)
  dev.off()
}

merge_nano_prot <- function(nano, prot) {
  nano <- nano$all
  prot <- prot$all
  df <- merge(nano, prot, by="pathway", suffixes=c(".nano", ".prot"))
  
  df$intersect <- mapply(intersect, df$leadingEdge.Ensembl.nano, df$leadingEdge.Ensembl.prot, SIMPLIFY=FALSE)
  df$size.intersect <- mapply(length, df$intersect)
  df$overlap_coefficient <- sign(df$NES.nano) * sign(df$NES.prot) * df$size.intersect / pmin(mapply(length, df$leadingEdge.Ensembl.nano), mapply(length, df$leadingEdge.GeneName.prot))

  return(df)
}

merge_two_models <- function(x, y, on='Ensembl', suffixes=c('.x','.y')) {
  #on can be 'Ensembl' or 'GeneName'
  df <- merge(x$all,y$all, by='pathway', suffixes=suffixes)
  
  NES.x <- paste('NES',suffixes[1],sep='')
  NES.y <- paste('NES',suffixes[2],sep='')
  leadingEdgeEnsembl.x <- paste('leadingEdge.Ensembl', suffixes[1],sep='')
  leadingEdgeGeneName.x <- paste('leadingEdge.GeneName', suffixes[1],sep='')
  leadingEdgeEnsembl.y <- paste('leadingEdge.Ensembl',suffixes[2],sep='')
  leadingEdgeGeneName.y <- paste('leadingEdge.GeneName',suffixes[2],sep='')
  pval.x <- paste('pval',suffixes[1],sep='')
  pval.y <- paste('pval',suffixes[2],sep='')
  
  df$intersect.Ensembl <- mapply(intersect, df[[leadingEdgeEnsembl.x]], df[[leadingEdgeEnsembl.y]], SIMPLIFY=FALSE)
  df$intersect.GeneName <- mapply(intersect, df[[leadingEdgeGeneName.x]], df[[leadingEdgeGeneName.y]], SIMPLIFY=FALSE)
  if (on=='Ensembl') {
    intersect_col = df$intersect.Ensembl
    leadingEdge.x <- leadingEdgeEnsembl.x
    leadingEdge.y <- leadingEdgeEnsembl.y
    } else {
    intersect_col = df$intersect.GeneName
    leadingEdge.x <- leadingEdgeGeneName.x
    leadingEdge.y <- leadingEdgeGeneName.y
    }
  
  df$size.intersect <- mapply(length, intersect_col)
  df$overlap_coefficient <- df$size.intersect / pmin(mapply(length, df[[leadingEdge.x]]), mapply(length, df[[leadingEdge.y]]))

  return(df)
}

plot_merged_gsea <- function(df, gs, path, filename, x_label, y_label, x_name, y_name) { 
  
  ## --- Column discovery (two runs) ---
  pval_cols <- grep("^pval\\.", names(df), value = TRUE)
  suffixes  <- unique(sub("^pval\\.", "", pval_cols))
  stopifnot(length(suffixes) == 2)
  
  suf_x <- suffixes[1]
  suf_y <- suffixes[2]
  
  NES_x <- paste0("NES.",  suf_x)
  NES_y <- paste0("NES.",  suf_y)
  FDR_x <- paste0("padj.", suf_x)
  FDR_y <- paste0("padj.", suf_y)
  
  ## --- Labels & tiers ---
  both   <- "FDR < 0.05 in both"
  only_x <- paste0(x_name, " FDR < 0.05 Only")
  only_y <- paste0(y_name, " FDR < 0.05 Only")
  
  df$pathway <- wrap.it(df$pathway, 25)
  stopifnot(length(df$pathway) > 0)
  
  ## --- Indicator rows (any ind.* TRUE) ---
  ind_cols <- grep("^ind\\.", names(df), value = TRUE)
  ind_any  <- if (length(ind_cols)) rowSums(df[, ..ind_cols] == TRUE, na.rm = TRUE) > 0 else rep(FALSE, nrow(df))
  
  ## --- FDR handling (robust numeric) ---
  fdr_x   <- suppressWarnings(as.numeric(as.character(df[[FDR_x]])))
  fdr_y   <- suppressWarnings(as.numeric(as.character(df[[FDR_y]])))
  fdr_min <- pmin(fdr_x, fdr_y, na.rm = TRUE); fdr_min[is.infinite(fdr_min)] <- NA_real_
  
  sig_any <- (fdr_x < 0.05) | (fdr_y < 0.05)
  sig_any[is.na(sig_any)] <- FALSE
  
  ## --- Pick labels: up to 10 per quadrant among significant rows, ind_any first ---
  df$label <- NA_character_
  
  pick_labels <- function(idx, n = 10) {
    if (!length(idx)) return(integer(0))
    idx_ind <- idx[ind_any[idx]]
    chosen  <- head(idx_ind[order(fdr_min[idx_ind], na.last = TRUE)], n)
    if (length(chosen) < n) {
      idx_rest <- setdiff(idx, chosen)
      chosen <- c(chosen, head(idx_rest[order(fdr_min[idx_rest], na.last = TRUE)], n - length(chosen)))
    }
    chosen
  }
  
  nes_x <- df[[NES_x]]
  nes_y <- df[[NES_y]]
  
  quadrants <- list(
    pp = which(sig_any & nes_x > 0 & nes_y > 0),
    pn = which(sig_any & nes_x > 0 & nes_y < 0),
    np = which(sig_any & nes_x < 0 & nes_y > 0),
    nn = which(sig_any & nes_x < 0 & nes_y < 0)
  )
  
  for (idx in quadrants) {
    chosen <- pick_labels(idx, n = 10)
    df$label[chosen] <- df$pathway[chosen]
  }
  
  ## --- Label repel nudges (radial outward) ---
  df_lab <- df[!is.na(label), ]
  
  r <- sqrt(df_lab[[NES_x]]^2 + df_lab[[NES_y]]^2)
  r[r == 0] <- NA_real_
  
  nudge_scale     <- 0.4
  df_lab$nudge_x  <- (df_lab[[NES_x]] / r) * nudge_scale
  df_lab$nudge_y  <- (df_lab[[NES_y]] / r) * nudge_scale
  
  ## --- Significance tier for outline legend ---
  df$sig_tier <- ifelse(fdr_x < 0.05 & fdr_y < 0.05, both,
                        ifelse(fdr_x < 0.05 & fdr_y >= 0.05, only_x,
                               ifelse(fdr_x >= 0.05 & fdr_y < 0.05, only_y, NA)))
  
  df$sig_tier <- factor(df$sig_tier, levels = c(only_x, only_y, both))
  
  outline_cols <- setNames(c("blue", "red", "purple"), c(only_x, only_y, both))
  
  ## --- Plot ---
  p <- ggplot(df, aes(x = .data[[NES_x]], y = .data[[NES_y]])) +
    
    # not significant in either
    geom_point(
      data = df[fdr_x >= 0.05 & fdr_y >= 0.05, ],
      aes(size = size.intersect, fill = overlap_coefficient),
      shape = 21, stroke = 0, alpha = 0.5
    ) + # significant in x only
    geom_point(
      data = df[fdr_x < 0.05 & fdr_y >= 0.05, ],
      aes(size = size.intersect, fill = overlap_coefficient, color = sig_tier),
      shape = 21, stroke = 0.5, alpha = 0.75
    ) + # significant in y only
    geom_point(
      data = df[fdr_x >= 0.05 & fdr_y < 0.05, ],
      aes(size = size.intersect, fill = overlap_coefficient, color = sig_tier),
      shape = 21, stroke = 0.5, alpha = 0.75
    ) + # significant in both
    geom_point(
      data = df[fdr_x < 0.05 & fdr_y < 0.05, ],
      aes(size = size.intersect, fill = overlap_coefficient, color = sig_tier),
      shape = 21, stroke = .7, alpha = 0.75
    ) +
    
    scale_color_manual(
      name = "BH-adjusted p (FDR)",
      values = outline_cols,
      drop = FALSE,
      na.translate = FALSE
    ) +
    guides(
      size = guide_legend(override.aes = list(color = "black", stroke = 0.8, alpha = 1))
    ) +
    
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    
    labs(
      x = paste0(x_label, " NES"),
      y = paste0(y_label, " NES"),
      title = paste0("GSEA results for ", x_name, " and ", y_name),
      size = "Intersection size",
      fill = "Overlap coefficient"
    ) +
    
    theme_light() +
    theme(
      plot.title = element_text(hjust = 0.5)#,
    #  plot.margin = margin(10, 30, 10, 30)
    ) +
    
    #coord_cartesian(clip = "off") +
    scale_x_continuous(expand = expansion(mult = 0.08)) +
    scale_y_continuous(expand = expansion(mult = 0.08)) +
    
    ggrepel::geom_label_repel(
      data = df_lab,
      aes(x = .data[[NES_x]], y = .data[[NES_y]], label = label),
      nudge_x = df_lab$nudge_x,
      nudge_y = df_lab$nudge_y,
      force_pull = 0,
      force = 10,
      box.padding = 0.5,
      point.padding = 0,
      min.segment.length = 0.1,
      max.time = 60,
      size = 1,
      alpha = 0.6,
      color = "black",
      fill = "white"
    )
  
    png(filename=paste(path, "/", filename, " merged GSEA plot", ".png", sep=""), width=2500, height=1700, res=300)
    print(p)
    dev.off()
    ggsave(paste(path, "/", filename, " merged GSEA plot", ".svg", sep=""), p, width = 8.33, height = 5.667)
    
}

if (!exists("mirna_gs", envir = globalenv())) {
msigdb_c3_mirna <- msigdbr(species="human",
                           collection="C3",
                           subcollection="MIR:MIRDB")
mirna_gs <- convert_to_genesets(msigdb_c3_mirna, "ensembl")
mirna_gs_gsymbol <- convert_to_genesets(msigdb_c3_mirna, "gene symbol")
}

if (!exists("tft_gs", envir = globalenv())) {
msigdb_c3_gtrd <- msigdbr(species="human",
                          collection="C3",
                          subcollection="GTRD")
msigdb_c3_legacy <- msigdbr(species="human",
                            collection="C3",
                            subcollection="TFT_LEGACY")
tft_gs <- combine_and_filter_genesets(msigdb_c3_gtrd, msigdb_c3_legacy) %>% convert_to_genesets("ensembl")
tft_gs_gsymbol <- combine_and_filter_genesets(msigdb_c3_gtrd, msigdb_c3_legacy) %>% convert_to_genesets("gene symbol")
}

if (!exists("h_gs", envir = globalenv())) {
msigdb_h <- msigdbr(species="human",
                    collection="H")
h_gs <- convert_to_genesets(msigdb_h, "ensembl")
h_gs_gsymbol <- convert_to_genesets(msigdb_h, "gene symbol")
}

if (!exists("reactome_gs", envir = globalenv())) {
reactome <- msigdbr(species="human",
                    collection="C2",
                    subcollection="CP:REACTOME")
reactome_gs <- convert_to_genesets(reactome, "ensembl")
reactome_gs_gsymbol <- convert_to_genesets(reactome, "gene symbol")
}

if (!exists("kegg_gs", envir = globalenv())) {
kegg_legacy <- msigdbr(species="human",
                      collection="C2",
                      subcollection="CP:KEGG_LEGACY")
kegg_gs <- convert_to_genesets(kegg_legacy, "ensembl")
kegg_gs_gsymbol <- convert_to_genesets(kegg_legacy, "gene symbol")
}

if (!exists("gocc_gs", envir = globalenv())) {
msigdb_go_cc <- msigdbr(species="human",
                     collection="C5",
                     subcollection="GO:CC")
gocc_gs <- convert_to_genesets(msigdb_go_cc, "ensembl")
gocc_gs_gsymbol <- convert_to_genesets(msigdb_go_cc, "gene symbol")
}

if (!exists("gobp_gs", envir = globalenv())) {
msigdb_go_bp <- msigdbr(species="human",
                     collection="C5",
                     subcollection="GO:BP")
gobp_gs <- convert_to_genesets(msigdb_go_bp, "ensembl")
gobp_gs_gsymbol <- convert_to_genesets(msigdb_go_bp, "gene symbol")
}

if (!exists("gomf_gs", envir = globalenv())) {
msigdb_go_mf <- msigdbr(species="human",
                     collection="C5",
                     subcollection="GO:MF")
gomf_gs <- convert_to_genesets(msigdb_go_mf, "ensembl")
gomf_gs_gsymbol <- convert_to_genesets(msigdb_go_mf, "gene symbol")
}

