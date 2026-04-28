library(tidyverse)
library(survival)
library(survminer)
library(readxl)
library(ggplot2)
library(cowplot)


outdir <- "/path/to//survival_analysis_output/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

target_pathways <- c(
  "GO_AMINOGLYCAN_BIOSYNTHETIC_PROCESS",
  "REACTOME_METABOLISM_OF_CARBOHYDRATES",
  "GO_GLYCOPROTEIN_METABOLIC_PROCESS",
  "REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM",
  "REACTOME_A_TETRASACCHARIDE_LINKER_SEQUENCE_IS_REQUIRED_FOR_GAG_SYNTHESIS"
)


ssgsea_file <- "/path/to//ssgsea_scores_all_glyco_pathways.csv"
ssgsea_scores <- read.csv(ssgsea_file, row.names = 1, check.names = FALSE)
rownames(ssgsea_scores) <- trimws(rownames(ssgsea_scores))



metadata_file <- "/path/to//cluster_assignments_with_cancertype.csv"
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
metadata$Sample <- trimws(gsub("\\.", "-", metadata$Sample))

clinical_file <- "/Users/arah/Downloads/clinic_pan_cancer_metada.xlsx"
clinical_data_all <- read_xlsx(clinical_file)
colnames(clinical_data_all) <- tolower(colnames(clinical_data_all))

required_cols <- c("bcr_patient_barcode", "type", "os", "os.time")
missing <- setdiff(required_cols, colnames(clinical_data_all))

clinical_data_all <- clinical_data_all %>%
  mutate(
    sample = trimws(gsub("\\.", "-", tolower(bcr_patient_barcode))),
    cancertype = trimws(tolower(type)),
    os = as.integer(os),
    os.time = as.numeric(os.time)
  ) %>%
  filter(!is.na(os.time), os.time > 0, !is.na(os)) %>%
  distinct(sample, .keep_all = TRUE)


prepare_pathway_df <- function(pathway, ssgsea_mat, meta, clinical) {
  if (!(pathway %in% rownames(ssgsea_mat))) {
    warning(paste("Pathway NOT FOUND:", pathway))
    return(data.frame())
  }
  
  scores <- ssgsea_mat[pathway, ]
  df <- data.frame(
    sample = colnames(ssgsea_mat),
    score = as.numeric(scores),
    stringsAsFactors = FALSE
  )

  df$sample_clean <- tolower(trimws(gsub("\\.", "-", df$sample)))
  meta$sample_clean <- tolower(trimws(gsub("\\.", "-", meta$Sample)))
  
  meta_unique <- meta %>% distinct(sample_clean, .keep_all = TRUE)
  
  df <- df %>%
    left_join(meta_unique[, c("sample_clean", "CancerType")], by = "sample_clean") %>%
    rename(cancertype = CancerType) %>%
    left_join(clinical[, c("sample", "os", "os.time")], by = c("sample_clean" = "sample")) %>%
    filter(!is.na(os), !is.na(os.time), !is.na(cancertype), !is.na(score)) %>%
    mutate(cancertype = tolower(trimws(cancertype)))
  
  df$os.time <- as.numeric(df$os.time)
  cat("  ✓", pathway, "→ n =", nrow(df), "\n")
  return(df)
}

pathway_data_list <- lapply(target_pathways, function(p) {
  prepare_pathway_df(p, ssgsea_scores, metadata, clinical_data_all)
})
names(pathway_data_list) <- target_pathways

for (i in seq_along(pathway_data_list)) {
  pathway_data_list[[i]]$pathway <- names(pathway_data_list)[i]
}


get_surv_categorized <- function(df) {
  out_list <- list()
  pathway_name <- unique(df$pathway)[1]
  
  for (ct in unique(df$cancertype)) {
    tmp <- df %>%
      filter(cancertype == ct) %>%
      mutate(os = as.integer(os), os.time = as.numeric(os.time)) %>%
      filter(!is.na(os.time), os.time > 0)
    
    if (nrow(tmp) < 30) {
      out_list[[ct]] <- data.frame()
      next
    }

    cut <- tryCatch({
      surv_cutpoint(tmp, time = "os.time", event = "os", variables = "score", minprop = 0.15)
    }, error = function(e) NULL)
    
    if (!is.null(cut)) {
      classed <- surv_categorize(cut)
      tmp$group_binary <- classed$score
      tmp$cutoff_method <- "surv_cutpoint"
    } else {

      median_cut <- median(tmp$score, na.rm = TRUE)
      tmp$group_binary <- ifelse(tmp$score > median_cut, "High", "Low")
      tmp$cutoff_method <- "median_split"
    }
    
    out_list[[ct]] <- tmp
  }
  
  bind_rows(out_list)
}

pathway_data_grouped <- lapply(pathway_data_list, get_surv_categorized)

perform_cox_analysis <- function(pathway_df, pathway_name) {
  cox_results <- list()
  
  for (cancer in unique(pathway_df$cancertype)) {
    cancer_data <- pathway_df %>%
      filter(cancertype == cancer, !is.na(group_binary))
    
    if (nrow(cancer_data) < 20 | length(unique(cancer_data$group_binary)) < 2) {
      next
    }
    
    tryCatch({
      cox_fit <- coxph(Surv(os.time, os) ~ group_binary, data = cancer_data)
      cox_summary <- summary(cox_fit)
      
      cox_results[[cancer]] <- data.frame(
        Pathway = pathway_name,
        CancerType = cancer,
        HR = cox_summary$conf.int[1, 1],
        HR_lower = cox_summary$conf.int[1, 3],
        HR_upper = cox_summary$conf.int[1, 4],
        pvalue = cox_summary$coefficients[1, 5],
        n = nrow(cancer_data),
        Cutoff_Method = unique(cancer_data$cutoff_method)
      )
      
      cat("Cox ✓:", pathway_name, "×", cancer, 
          "HR =", round(cox_summary$conf.int[1, 1], 2),
          "p =", format.pval(cox_summary$coefficients[1, 5], digits = 2), "\n")
    }, error = function(e) {
      cat("Cox ✗:", pathway_name, "×", cancer, "-", e$message, "\n")
    })
  }
  
  bind_rows(cox_results)
}

cox_results_list <- Map(perform_cox_analysis, pathway_data_grouped, target_pathways)
all_cox_results <- bind_rows(cox_results_list)



for (i in seq_along(cox_results_list)) {
  pathway_name <- names(cox_results_list)[i]
  cox_df <- cox_results_list[[i]]
  
  if (nrow(cox_df) == 0) next
  
  cox_df <- cox_df %>%
    mutate(sig = ifelse(pvalue < 0.05, "***", ifelse(pvalue < 0.1, "*", "")))
  
  p <- ggplot(cox_df, aes(x = reorder(CancerType, HR), y = HR)) +
    geom_point(aes(color = ifelse(pvalue < 0.05, "Significant", "NS")), size = 5) +
    geom_errorbar(aes(ymin = HR_lower, ymax = HR_upper, color = ifelse(pvalue < 0.05, "Significant", "NS")),
                  width = 0.3, linewidth = 1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "darkred", linewidth = 1) +
    coord_flip() +
    scale_color_manual(name = "p-value", values = c("Significant" = "#E41A1C", "NS" = "#999999"),
                       labels = c("Significant" = "p < 0.05", "NS" = "p ≥ 0.05")) +
    scale_y_log10(breaks = c(0.5, 1, 2, 4)) +
    theme_bw(base_size = 14) +
    labs(
      title = gsub("_", " ", pathway_name),
      subtitle = "Cox Proportional Hazards - High vs Low Pathway Activity",
      x = "Cancer Type",
      y = "Hazard Ratio (95% CI, log scale)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  filename <- paste0(outdir, "cox_forest_", tolower(gsub(" ", "_", pathway_name)), ".tiff")
  ggsave(filename, p, device = "tiff", width = 10, height = 7, dpi = 600, compression = "lzw")
  cat("  ✓ Saved:", filename, "\n")
}


library(stringr)

beta_heatmap_data <- all_cox_results %>%
  filter(!is.na(HR)) %>%
  mutate(
    beta = log(HR),
    sig_label = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01  ~ "**",
      pvalue < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Pathway_label = str_replace_all(Pathway, "_", " ")
  )


beta_heatmap_data$Pathway_label <- ifelse(
  beta_heatmap_data$Pathway ==
    "REACTOME_A_TETRASACCHARIDE_LINKER_SEQUENCE_IS_REQUIRED_FOR_GAG_SYNTHESIS",
  "REACTOME A TETRASACCHARIDE LINKER\nSEQUENCE IS REQUIRED FOR GAG SYNTHESIS",
  beta_heatmap_data$Pathway_label
)

beta_heatmap_data$Pathway_label <- factor(
  beta_heatmap_data$Pathway_label,
  levels = rev(unique(beta_heatmap_data$Pathway_label))
)

p_heatmap <- ggplot(beta_heatmap_data,
                    aes(x = CancerType, y = Pathway_label, fill = beta)) +
  geom_tile(color = "grey90", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1), oob = scales::squish,
    name = "Log(HR)\n(Beta)"
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Cox Regression - Pathway Activity and Survival",
    x = "Cancer Type", y = "Pathway"
  ) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    panel.grid  = element_blank(),
    plot.margin = margin(5, 5, 5, 5, "pt")
  )

filename <- paste0("/path/to//survival_analysis_output/beta_coefficient_heatmap.tiff")
ggsave(filename, p_heatmap, device = "tiff", width = 12, height = 8, dpi = 600, compression = "lzw")


for (i in seq_along(pathway_data_grouped)) {
  pathway_name <- names(pathway_data_grouped)[i]
  pathway_df <- pathway_data_grouped[[i]]
  
  if (nrow(pathway_df) == 0) next
  
  p_violin <- ggplot(pathway_df, aes(x = cancertype, y = score, fill = cancertype)) +
    geom_violin(alpha = 0.7, trim = TRUE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 2) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw(base_size = 13) +
    labs(
      title = gsub("_", " ", pathway_name),
      subtitle = "Score Distribution across Cancer Types",
      x = "Cancer Type",
      y = "ssGSEA Enrichment Score"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  filename <- paste0(outdir, "score_distribution_", tolower(gsub(" ", "_", pathway_name)), ".tiff")
  ggsave(filename, p_violin, device = "tiff", width = 10, height = 7, dpi = 600, compression = "lzw")
  cat("  ✓ Saved:", filename, "\n")
}



for (i in seq_along(pathway_data_grouped)) {
  pathway_name <- names(pathway_data_grouped)[i]
  pathway_df <- pathway_data_grouped[[i]]
  
  for (cancer in unique(pathway_df$cancertype)) {
    plot_data <- pathway_df %>%
      filter(cancertype == cancer, !is.na(group_binary))
    
    if (nrow(plot_data) < 20 | length(unique(plot_data$group_binary)) < 2) {
      next
    }
    
    plot_data$os.time <- as.numeric(plot_data$os.time)
    fit <- survfit(Surv(os.time, os) ~ group_binary, data = plot_data)
    
    p_km <- ggsurvplot(
      fit, data = plot_data,
      pval = TRUE, pval.size = 5,
      conf.int = TRUE,
      risk.table = TRUE, risk.table.height = 0.3,
      legend.title = "Pathway Activity",
      legend.labs = c("Low", "High"),
      palette = c("Low" = "#2166AC", "High" = "#B2182B"),
      xlab = "Time (days)",
      ylab = "Overall Survival Probability",
      title = paste(gsub("_", " ", pathway_name), "-", toupper(cancer)),
      ggtheme = theme_bw(base_size = 13)
    )
    
    filename <- paste0(outdir, "kaplan_meier_", tolower(gsub(" ", "_", pathway_name)), "_", tolower(cancer), ".tiff")
    ggsave(filename, p_km$plot, device = "tiff", width = 10, height = 8, dpi = 600, compression = "lzw")
    cat("  ✓ Saved:", filename, "\n")
  }
}

summary_table <- all_cox_results %>%
  arrange(Pathway, pvalue) %>%
  mutate(
    HR_CI = paste0(format(round(HR, 2), nsmall = 2), 
                   " (", format(round(HR_lower, 2), nsmall = 2), 
                   "-", format(round(HR_upper, 2), nsmall = 2), ")"),
    pvalue_formatted = format.pval(pvalue, digits = 3),
    Significant = ifelse(pvalue < 0.05, "Yes", "No")
  ) %>%
  select(Pathway, CancerType, HR_CI, pvalue_formatted, Significant, n, Cutoff_Method)

write.csv(summary_table, paste0(outdir, "pathway_cox_summary_table.csv"), row.names = FALSE)



