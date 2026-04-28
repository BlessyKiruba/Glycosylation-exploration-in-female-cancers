
library(GSVA)
library(GSEABase)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(viridis)
library(readxl)

expr_file <- "/path/to/common_genesfromboth_iqr_mad.csv"
gmt_file  <- "/path/to/20231120_Glyco_gene_sets.gmt"
meta_file <- "/path/to/cluster_assignments_with_cancerType.csv"

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_file)

expr <- expr[, colnames(expr) %in% meta$Sample]
meta <- meta[match(colnames(expr), meta$Sample), ]

if (max(expr, na.rm = TRUE) > 100) {
  expr <- log2(expr + 1)
}

geneSets <- getGmt(gmt_file)

ssgsea_scores <- gsva(
  as.matrix(expr),
  geneSets,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = TRUE,
  min.sz = 10,
  max.sz = 5000,
  parallel.sz = 1,
  verbose = TRUE
)

ssgsea_scores_df <- as.data.frame(ssgsea_scores)
write.csv(
  ssgsea_scores_df,
  "/path/to/ssgsea_scores_all_glyco_pathways.csv",
  quote = FALSE
)

ssgsea_file <- "/path/to/ssgsea_scores_all_glyco_pathways.csv"
metadata_file <- "/path/to/cluster_assignments_with_cancertype.csv"

ssgsea <- read.csv(ssgsea_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file)

ssgsea_long <- ssgsea %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(
    cols = -Pathway,
    names_to = "Sample",
    values_to = "Score"
  )

merged_df <- ssgsea_long %>%
  left_join(metadata, by = "Sample") %>%
  filter(!is.na(CancerType))

summary_df <- merged_df %>%
  group_by(Pathway, CancerType) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop")

top30_df <- summary_df %>%
  group_by(CancerType) %>%
  slice_max(order_by = mean_score, n = 30) %>%
  ungroup()

pathway_counts <- top30_df %>%
  group_by(Pathway) %>%
  summarise(n_cancers = n_distinct(CancerType), .groups = "drop") %>%
  arrange(desc(n_cancers))

common_pathways <- pathway_counts %>%
  filter(n_cancers == 4)

selected_pathways <- c(
  "GO_AMINOGLYCAN_BIOSYNTHETIC_PROCESS",
  "GO_GLYCOPROTEIN_METABOLIC_PROCESS",
  "REACTOME_A_TETRASACCHARIDE_LINKER_SEQUENCE_IS_REQUIRED_FOR_GAG_SYNTHESIS",
  "REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM",
  "REACTOME_METABOLISM_OF_CARBOHYDRATES"
)

top30_df_filtered <- top30_df %>%
  filter(Pathway %in% selected_pathways)

bubble_plot_selected <- ggplot(
  top30_df_filtered,
  aes(
    x = CancerType,
    y = reorder(
      stringr::str_wrap(gsub("_", " ", Pathway), width = 32),
      mean_score
    )
  )
) +
  geom_point(
    aes(size = mean_score, fill = mean_score),
    shape = 21,
    color = "black",
    stroke = 0.4,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(
    option = "magma",
    begin = 0.15,
    end = 0.95
  ) +
  scale_size_continuous(
    range = c(4, 12)
  ) +
  scale_x_discrete(
    expand = expansion(add = 0.6)
  ) +
  scale_y_discrete(
    expand = expansion(add = 0.4)
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 16, 12, 12),
    legend.position = "right"
  ) +
  labs(
    title = "Top Glycosylation Pathways Across Cancer Types",
    subtitle = "Bubble size and color represent mean ssGSEA scores",
    x = "Cancer Type",
    y = "Pathway",
    fill = "Mean ssGSEA\nScore",
    size = "Mean ssGSEA\nScore"
  )

print(bubble_plot_selected)

ggsave(
  "/path/to/bubbleplot_top30_ssgsea_cancertype.tiff",
  bubble_plot_selected,
  width = 15,
  height = 7,
  dpi = 600
)

bubble_plot_top30 <- ggplot(
  top30_df,
  aes(
    x = CancerType,
    y = reorder(Pathway, mean_score)
  )
) +
  geom_point(
    aes(size = mean_score, color = mean_score),
    alpha = 0.8
  ) +
  scale_color_viridis_c(option = "magma") +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank()
  ) +
  labs(
    title = "Top 30 ssGSEA Pathways per Cancer Type",
    x = "Cancer Type",
    y = "Pathway",
    color = "Mean Score",
    size = "Mean Score"
  )

print(bubble_plot_top30)

ggsave(
  "/path/to/bubbleplot_top30_ssgsea_cancertype.png",
  bubble_plot_top30,
  width = 15,
  height = 10,
  dpi = 500
)
