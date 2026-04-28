library(dplyr)
library(ggplot2)
library(ggpubr)

df <- read.csv(
  "/path/to//Top_genes_Tumor_vs_Normal_pan_cancer.csv",
  stringsAsFactors = FALSE
)

df$CancerType <- factor(df$CancerType, levels = c("BRCA", "CC", "OV", "UCEC"))
df$SampleType <- factor(df$SampleType, levels = c("Normal", "Tumor"))
df$Gene       <- factor(df$Gene, levels = c("CP", "DAPL1", "MUC20", "SLC16A3", "TFF3"))

p <- ggplot(df, aes(x = SampleType, y = Expression_plot, fill = SampleType)) +
  geom_boxplot(width = 0.6, outlier.shape = 16, outlier.size = 0.8, alpha = 0.85,
               linewidth = 0.4) +
  stat_compare_means(
    method       = "wilcox.test",
    label        = "p.format",
    size         = 3.2,
    label.x.npc  = "left",
    label.y.npc  = "top"
  ) +
  facet_grid(CancerType ~ Gene, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#F4A582", "Tumor" = "#92C5DE")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    title = "Top genes: Tumor vs Normal across four cancers",
    x     = "Sample type",
    y     = "log2(count + 1)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position       = "none",
    strip.background      = element_rect(fill = "grey92", colour = "grey50", linewidth = 0.4),
    strip.text            = element_text(face = "bold", size = 11),
    axis.text.x           = element_text(angle = 45, hjust = 1, size = 10),
    axis.title            = element_text(size = 11),
    plot.title            = element_text(face = "bold", size = 13, hjust = 0),
    panel.grid.minor      = element_blank(),
    panel.border          = element_rect(colour = "grey70", linewidth = 0.4)
  )

ggsave(
  filename    = "/path/to//Top_genes_Tumor_vs_Normal_four_cancers_final.tiff",
  plot        = p,
  width       = 14,
  height      = 12,
  dpi         = 600,
  compression = "lzw",
  units       = "in"
)
