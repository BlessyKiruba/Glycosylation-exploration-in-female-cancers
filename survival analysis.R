
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(purrr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
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


data_merged <- inner_join(expr, clinic, by = "sample_id")


data_merged <- data_merged %>%
  arrange(sample_id) %>%
  distinct(sample_id, .keep_all = TRUE)


stopifnot(!any(duplicated(data_merged$sample_id)))
stopifnot(!any(is.na(data_merged$sample_id)))
stopifnot(all(c("CancerType", "time", "status") %in% colnames(data_merged)))

print(table(data_merged$CancerType))


run_survival_gene <- function(df, gene) {
  
  df <- df %>% filter(!is.na(.data[[gene]]))
  if (nrow(df) < 30) return(NULL)
  
  cut <- tryCatch(
    surv_cutpoint(
      df,
      time = "time",
      event = "status",
      variables = gene,
      minprop = 0.2
    ),
    error = function(e) NULL
  )
  
  if (is.null(cut)) return(NULL)
  
  cat_df <- surv_categorize(cut)
  
  if (min(table(cat_df[[gene]])) < 10) return(NULL)
  
  fit <- coxph(
    as.formula(paste0("Surv(time, status) ~ ", gene)),
    data = cat_df
  )
  
  s <- summary(fit)
  
  tibble(
    Gene = gene,
    HR = s$coefficients[, "exp(coef)"],
    pvalue = s$coefficients[, "Pr(>|z|)"],
    Direction = ifelse(HR < 1, "Higher Survival", "Lower Survival")
  )
}


gene_list <- setdiff(
  colnames(data_merged),
  c("sample_id", "CancerType", "time", "status")
)


plot_df <- surv_results %>%
  mutate(
    SurvivalGroup = ifelse(HR < 1, "Higher survival", "Lower survival"),
    negLogP = -log10(pvalue),
    CancerType = factor(CancerType, levels = c("BRCA", "CC", "OV", "UCEC"))
  )
p_high <- ggplot(
  plot_df %>% filter(SurvivalGroup == "Higher survival"),
  aes(
    x = CancerType,
    y = reorder(Gene, logHR),
    size = negLogP,
    fill = logHR
  )
) +
  geom_point(shape = 21, color = "black", stroke = 0.4) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "log2(HR)"
  ) +
  scale_size(range = c(2.5, 9), name = expression(-log[10](italic(p)))) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Higher survival–associated genes",
    x = "Cancer type",
    y = "Gene"
  )
p_low <- ggplot(
  plot_df %>% filter(SurvivalGroup == "Lower survival"),
  aes(
    x = CancerType,
    y = reorder(Gene, logHR),
    size = negLogP,
    fill = logHR
  )
) +
  geom_point(shape = 21, color = "black", stroke = 0.4) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "log2(HR)"
  ) +
  scale_size(range = c(2.5, 9), name = expression(-log[10](italic(p)))) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Lower survival–associated genes",
    x = "Cancer type",
    y = ""
  )
final_plot <- p_high | p_low
print(final_plot)


ggsave(
  "/path/to/glycogene_survival_bubbleplot_in_pancancer.tiff",
  final_plot,
  width = 14,
  height = 12,
  dpi = 600,
)

write.csv(
  data_merged,
  "/path/to/compiled_63gene_expression_survival.csv",
  row.names = FALSE
)

