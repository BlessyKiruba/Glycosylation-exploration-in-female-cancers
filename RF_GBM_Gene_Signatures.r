
rm(list = ls())
gc()

setwd("/path/to/")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  survival, survminer, randomForestSRC, gbm, caret, pROC,
  tidyverse, ggplot2, survcomp
)



data <- read.csv("/path/to/survival_expression_ml_compiled.csv", stringsAsFactors = FALSE)

cat(sprintf("Loaded: %d samples x %d columns\n", nrow(data), ncol(data)))

exclude_cols <- c("sample_id", "time", "status")
gene_cols <- setdiff(colnames(data), exclude_cols)

cat(sprintf("Number of gene features: %d\n", length(gene_cols)))
cat(sprintf("Events: %d/%d (%.1f%%)\n\n",
            sum(data$status), nrow(data),
            mean(data$status) * 100))

X <- as.matrix(data[, gene_cols])
time <- data$time
event <- data$status

surv_obj <- Surv(time, event)
X_scaled <- scale(X)


set.seed(42)

train_val_idx <- sample(1:nrow(X_scaled), size = 0.8 * nrow(X_scaled))
test_idx <- setdiff(1:nrow(X_scaled), train_val_idx)

train_idx <- sample(train_val_idx, size = 0.75 * length(train_val_idx))
val_idx <- setdiff(train_val_idx, train_idx)

X_train <- X_scaled[train_idx, ]
X_val   <- X_scaled[val_idx, ]
X_test  <- X_scaled[test_idx, ]

time_train <- time[train_idx]
time_val   <- time[val_idx]
time_test  <- time[test_idx]

event_train <- event[train_idx]
event_val   <- event[val_idx]
event_test  <- event[test_idx]

surv_train <- Surv(time_train, event_train)
surv_val   <- Surv(time_val, event_val)
surv_test  <- Surv(time_test, event_test)

cat(sprintf("Train: %d | Val: %d | Test: %d\n\n",
            length(train_idx), length(val_idx), length(test_idx)))


results <- data.frame(
  Model = character(),
  Val_C_Index = numeric(),
  Test_C_Index = numeric(),
  Status = character(),
  stringsAsFactors = FALSE
)

rsf_rank <- read.csv("feature_selection_rsf_importance.csv", stringsAsFactors = FALSE)

rsf_topk_genes <- head(rsf_rank$Gene, 26)
gbm_topk_genes <- head(rsf_rank$Gene, 43)


tryCatch({
  rsf_genes <- intersect(rsf_topk_genes, gene_cols)
  X_train_rsf <- X_train[, rsf_genes, drop = FALSE]
  X_val_rsf   <- X_val[,   rsf_genes, drop = FALSE]
  X_test_rsf  <- X_test[,  rsf_genes, drop = FALSE]
  
  train_df_rsf <- data.frame(X_train_rsf, time = time_train, event = event_train)
  
  rsf_model <- rfsrc(
    Surv(time, event) ~ .,
    data = train_df_rsf,
    ntree = 300,
    nodesize = 5,
    importance = "permute"
  )
  
  val_df_rsf  <- data.frame(X_val_rsf)
  test_df_rsf <- data.frame(X_test_rsf)
  
  rsf_pred_val  <- predict(rsf_model, newdata = val_df_rsf)$predicted
  rsf_pred_test <- predict(rsf_model, newdata = test_df_rsf)$predicted
  
  rsf_cindex_val <- survcomp::concordance.index(
    rsf_pred_val, surv.time = time_val, surv.event = event_val
  )$c.index
  
  rsf_cindex_test <- survcomp::concordance.index(
    rsf_pred_test, surv.time = time_test, surv.event = event_test
  )$c.index
  
  results <- rbind(results, data.frame(
    Model = "RSF",
    Val_C_Index = rsf_cindex_val,
    Test_C_Index = rsf_cindex_test,
    Status = "Success"
  ))
  
  cat(sprintf("  Val C-Index: %.4f | Test C-Index: %.4f\n",
              rsf_cindex_val, rsf_cindex_test))
  
  rsf_importance <- rsf_model$importance
  rsf_importance_df <- data.frame(
    Gene = names(rsf_importance),
    Importance = as.numeric(rsf_importance),
    stringsAsFactors = FALSE
  )
  
  rsf_importance_df <- rsf_importance_df[order(-rsf_importance_df$Importance), ]
  rsf_top10 <- head(rsf_importance_df, 10)
  
  cat("  Top 10 genes by RSF importance (within top-26 feature set):\n")
  print(rsf_top10)
  cat("\n")
  
  write.csv(rsf_importance_df, "rsf_feature_importance_top26_genes.csv", row.names = FALSE)
  write.csv(rsf_top10, "rsf_top10_gene_signatures_top26.csv", row.names = FALSE)

  
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n\n", e$message))
})


tryCatch({
  gbm_genes <- intersect(gbm_topk_genes, gene_cols)
  X_train_gbm <- X_train[, gbm_genes, drop = FALSE]
  X_val_gbm   <- X_val[,   gbm_genes, drop = FALSE]
  X_test_gbm  <- X_test[,  gbm_genes, drop = FALSE]
  
  train_df_gbm <- data.frame(X_train_gbm, time = time_train, event = event_train)
  
  gbm_model <- gbm(
    Surv(time, event) ~ .,
    data = train_df_gbm,
    distribution = "coxph",
    n.trees = 300,
    interaction.depth = 4,
    shrinkage = 0.001,
    n.minobsinnode = 10
  )
  
  val_df_gbm  <- data.frame(X_val_gbm)
  test_df_gbm <- data.frame(X_test_gbm)
  
  gbm_pred_val <- predict(gbm_model, newdata = val_df_gbm, n.trees = 300, type = "link")
  gbm_pred_test <- predict(gbm_model, newdata = test_df_gbm, n.trees = 300, type = "link")
  
  gbm_cindex_val <- survcomp::concordance.index(
    gbm_pred_val, surv.time = time_val, surv.event = event_val
  )$c.index
  
  gbm_cindex_test <- survcomp::concordance.index(
    gbm_pred_test, surv.time = time_test, surv.event = event_test
  )$c.index
  
  results <- rbind(results, data.frame(
    Model = "GBM",
    Val_C_Index = gbm_cindex_val,
    Test_C_Index = gbm_cindex_test,
    Status = "Success"
  ))
  
  cat(sprintf("  Val C-Index: %.4f | Test C-Index: %.4f\n",
              gbm_cindex_val, gbm_cindex_test))
  
  gbm_importance <- summary(gbm_model, plotit = FALSE)
  gbm_importance_df <- data.frame(
    Gene = gbm_importance$var,
    Importance = gbm_importance$rel.inf,
    stringsAsFactors = FALSE
  )
  
  gbm_importance_df <- gbm_importance_df[order(-gbm_importance_df$Importance), ]
  gbm_top10 <- head(gbm_importance_df, 20)
  
  cat("  Top 10 genes by GBM importance (within top-43 feature set):\n")
  print(gbm_top10)
  cat("\n")
  
  write.csv(gbm_importance_df, "gbm_feature_importance_top43_genes.csv", row.names = FALSE)
  write.csv(gbm_top10, "gbm_top10_gene_signatures_top43.csv", row.names = FALSE)

}, error = function(e) {
  cat(sprintf("  ERROR: %s\n\n", e$message))
})


rsf_top10_genes <- rsf_top10$Gene
gbm_top10_genes <- gbm_top10$Gene

cat(sprintf("RSF Top 10 genes (top-26 set): %s\n", paste(rsf_top10_genes, collapse = ", ")))
cat(sprintf("GBM Top 10 genes (top-43 set): %s\n\n", paste(gbm_top10_genes, collapse = ", ")))

common_genes <- intersect(rsf_top10_genes, gbm_top10_genes)

cat(sprintf("Number of common genes (consensus signatures): %d\n\n", length(common_genes)))

print(common_genes)
cat("\n")

if (length(common_genes) > 0) {
  consensus_df <- data.frame(
    Gene = common_genes,
    RSF_Rank = match(common_genes, rsf_top10_genes),
    RSF_Importance = rsf_top10$Importance[match(common_genes, rsf_top10_genes)],
    GBM_Rank = match(common_genes, gbm_top10_genes),
    GBM_Importance = gbm_top10$Importance[match(common_genes, gbm_top10_genes)],
    stringsAsFactors = FALSE
  )
  
  consensus_df$Avg_Rank <- (consensus_df$RSF_Rank + consensus_df$GBM_Rank) / 2
  consensus_df <- consensus_df[order(consensus_df$Avg_Rank), ]
  
  cat("Consensus signatures with rankings:\n")
  print(consensus_df)
  cat("\n")
  
  write.csv(consensus_df, "consensus_gene_signatures_rsf_gbm_topk.csv", row.names = FALSE)

}



results <- results[order(-results$Test_C_Index), ]
rownames(results) <- NULL

print(results)

best_model <- results[1, ]
cat("\n")
cat(sprintf("BEST MODEL: %s (Test C-Index: %.4f)\n",
            best_model$Model, best_model$Test_C_Index))



write.csv(results, "rf_gbm_topk_model_results.csv", row.names = FALSE)

rsf_top10_plot <- rsf_top10 %>%
  mutate(Gene = factor(Gene, levels = rev(Gene))) %>%
  ggplot(aes(x = Importance, y = Gene, fill = Importance)) +
  geom_col(color = "black", size = 0.8) +
  scale_fill_gradient(low = "#38bdf8", high = "#2ecc71") +
  labs(
    title = "Random Survival Forest - Top 10 Gene Signatures",
    subtitle = "Feature Importance for Pan-Cancer Survival Prediction",
    x = "Importance Score",
    y = "Gene"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

ggsave("rsf_top10_gene_signatures_plot.tiff", rsf_top10_plot, width = 10, height = 7, dpi = 600)
cat("Saved: rsf_top10_gene_signatures_plot.tiff\n")

# GBM Top 10 plot
gbm_top10_plot <- gbm_top10 %>%
  mutate(Gene = factor(Gene, levels = rev(Gene))) %>%
  ggplot(aes(x = Importance, y = Gene, fill = Importance)) +
  geom_col(color = "black", size = 0.8) +
  scale_fill_gradient(low = "#38bdf8", high = "#e74c3c") +
  labs(
    title = "Gradient Boosting Machine - Top 10 Gene Signatures",
    subtitle = "Feature Importance for Pan-Cancer Survival Prediction",
    x = "Importance Score (Relative Influence)",
    y = "Gene"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

ggsave("gbm_top10_gene_signatures_plot.tiff", gbm_top10_plot, width = 10, height = 7, dpi = 600)

comparison_plot <- data.frame(
  Gene = c(rsf_top10_genes, gbm_top10_genes),
  Model = c(rep("RSF", 10), rep("GBM", 10)),
  Rank = c(1:10, 1:10)
) %>%
  ggplot(aes(x = Model, y = Rank, fill = Gene)) +
  geom_tile(color = "black", size = 0.5) +
  scale_y_reverse() +
  labs(
    title = "RSF vs GBM - Top 10 Gene Signatures Comparison",
    subtitle = "Gene overlap (consensus signatures in intersection)",
    x = "Model",
    y = "Rank"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

ggsave("rsf_gbm_gene_comparison.tiff", comparison_plot, width = 10, height = 8, dpi = 600)
cat("Saved: rsf_gbm_gene_comparison.tiff\n")

