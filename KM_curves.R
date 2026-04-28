
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

expr <- read.csv(
  "/path/to/survival_expression_ml_compiled.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expr$sample_id <- gsub("\\.", "-", expr$sample_id)
expr$sample_id <- substr(expr$sample_id, 1, 12)

clinic <- read.csv(
  "/path/to/cluster_assignments_with_cancertype.csv",
  stringsAsFactors = FALSE
)

clinic$sample_id <- gsub("\\.", "-", clinic$sample_id)
clinic$sample_id <- substr(clinic$sample_id, 1, 12)

data_merged <- inner_join(expr, clinic, by = "sample_id") %>%
  arrange(sample_id) %>%
  distinct(sample_id, .keep_all = TRUE)

stopifnot(!any(duplicated(data_merged$sample_id)))
stopifnot(!any(is.na(data_merged$sample_id)))
stopifnot(all(c("CancerType", "time", "status") %in% colnames(data_merged)))

print(table(data_merged$CancerType))

genes_5 <- c("MUC20", "DAPL1", "SLC16A3", "TFF3", "CP")

median_split_df <- function(df, gene, min_n = 50, min_group_n = 10) {
  if (!gene %in% colnames(df)) return(NULL)
  
  d <- df %>%
    select(sample_id, time, status, all_of(gene)) %>%
    rename(expr = all_of(gene)) %>%
    filter(!is.na(expr), !is.na(time), !is.na(status))
  
  if (nrow(d) < min_n) return(NULL)
  
  med <- median(d$expr, na.rm = TRUE)
  if (!is.finite(med)) return(NULL)
  
  d <- d %>%
    mutate(expr_group = ifelse(expr > med, "High", "Low")) %>%
    mutate(expr_group = factor(expr_group, levels = c("Low", "High")))
  
  if (min(table(d$expr_group)) < min_group_n) return(NULL)
  
  d
}

#-------------------------------#
# Pan-cancer KM by gene
#-------------------------------#
plot_km_pan_gene_median <- function(df, gene) {
  d <- median_split_df(df, gene, min_n = 50, min_group_n = 10)
  if (is.null(d)) return(NULL)
  
  fit <- survfit(Surv(time, status) ~ expr_group, data = d)
  
  ggsurvplot(
    fit,
    data = d,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.28,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    legend.title = gene,
    legend.labs = c("Low (≤ median)", "High (> median)"),
    palette = c("#2166ac", "#b2182b"),
    title = paste0(gene, " (Pan-cancer, median split)"),
    xlab = "Time",
    ylab = "Overall survival probability",
    ggtheme = theme_bw(base_size = 13),
    risk.table.theme = theme_bw(base_size = 11),
    surv.median.line = "none"  
  )
}

km_pan_list <- lapply(genes_5, function(g) plot_km_pan_gene_median(data_merged, g))
names(km_pan_list) <- genes_5
km_pan_list <- km_pan_list[!sapply(km_pan_list, is.null)]

for (g in names(km_pan_list)) {
  surv_obj <- km_pan_list[[g]]
  
  combined <- arrange_ggsurvplots(
    list(surv_obj),
    ncol = 1,
    nrow = 1,
    print = FALSE
  )
  
  outfile <- paste0("/path/to/KM_pan_cancer_", g, "_median.tiff")
  
  ggsave(
    filename = outfile,
    plot = combined,
    width = 7,
    height = 8,
    dpi = 600,
    compression = "lzw"
  )
  
  message("Saved: ", outfile)
}

plot_km_cancer_gene_median <- function(df, gene, cancer) {
  d <- median_split_df(df, gene, min_n = 30, min_group_n = 10)
  if (is.null(d)) return(NULL)
  
  fit <- survfit(Surv(time, status) ~ expr_group, data = d)
  
  ggsurvplot(
    fit,
    data = d,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.28,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    legend.title = gene,
    legend.labs = c("Low (≤ median)", "High (> median)"),
    palette = c("#2166ac", "#b2182b"),
    title = paste0(gene, " (", cancer, ", median split)"),
    xlab = "Time",
    ylab = "Overall survival probability",
    ggtheme = theme_bw(base_size = 13),
    risk.table.theme = theme_bw(base_size = 11),
    surv.median.line = "none"
  )
}

outdir <- "/path/to/KM curves/median_split"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cancer_types <- sort(unique(data_merged$CancerType))

for (ct in cancer_types) {
  df_ct <- data_merged %>% filter(CancerType == ct)
  
  for (g in genes_5) {
    surv_obj <- plot_km_cancer_gene_median(df_ct, g, ct)
    if (is.null(surv_obj)) next
    
    combined <- arrange_ggsurvplots(
      list(surv_obj),
      ncol = 1,
      nrow = 1,
      print = FALSE
    )
    
    outfile <- file.path(outdir, paste0("KM_", ct, "_", g, "_median.tiff"))
    
    ggsave(
      filename = outfile,
      plot = combined,
      width = 7,
      height = 8,
      dpi = 600,
      compression = "lzw"
    )
    
    message("Saved: ", outfile)
  }
}
