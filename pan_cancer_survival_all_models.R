
rm(list = ls())
gc()

setwd("/path/to/")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  survival, glmnet, MASS, CoxBoost, randomForestSRC,
  plsRcox, superpc, gbm, tidyverse, ggplot2, survcomp
)


data1 <- read.csv("data_ridge_stepwise_filtered.csv")
data2 <- read.csv("data_lasso_elasticnet_filtered.csv")
data3 <- read.csv("data_coxboost_filtered.csv")
data4 <- read.csv("data_rsf_gbm_filtered.csv")
data5 <- read.csv("data_pca_plsrcox_superpc_filtered.csv")

rank1 <- read.csv("feature_selection_univariate_cox.csv")
rank2 <- read.csv("feature_selection_lasso_coefficients.csv")
rank3 <- read.csv("feature_selection_univariate_cox.csv")
rank4 <- read.csv("feature_selection_rsf_importance.csv")
rank5 <- read.csv("feature_selection_pca_loadings.csv")

time_raw  <- data1$time
event_raw <- data1$status

min_time <- min(time_raw, na.rm = TRUE)
eps <- 1e-6
shift <- if (min_time <= 0) abs(min_time) + eps else 0
time_fixed <- time_raw + shift

set.seed(42)
n <- length(time_fixed)
train_val_idx <- sample(1:n, size = 0.8 * n)
test_idx  <- setdiff(1:n, train_val_idx)
train_idx <- sample(train_val_idx, size = 0.75 * length(train_val_idx))
val_idx   <- setdiff(train_val_idx, train_idx)

time_train <- time_fixed[train_idx]
time_val   <- time_fixed[val_idx]
time_test  <- time_fixed[test_idx]

event_train <- event_raw[train_idx]
event_val   <- event_raw[val_idx]
event_test  <- event_raw[test_idx]

Surv_train <- Surv(time_train, event_train)

cat(sprintf("[SETUP] Common split: Train %d | Val %d | Test %d\n\n",
            length(train_idx), length(val_idx), length(test_idx)))


all_topk_results <- data.frame(
  Model = character(),
  Dataset = character(),
  k = integer(),
  Val_C = numeric(),
  Test_C = numeric(),
  stringsAsFactors = FALSE
)


fit_ridge_topk <- function(genes, dat) {
  genes <- intersect(genes, colnames(dat))
  if (length(genes) == 0) return(c(NA, NA))
  
  X_train <- scale(as.matrix(dat[train_idx, genes, drop = FALSE]))
  X_val   <- scale(as.matrix(dat[val_idx,   genes, drop = FALSE]))
  X_test  <- scale(as.matrix(dat[test_idx,  genes, drop = FALSE]))
  
  cv_ridge <- cv.glmnet(X_train, Surv_train, family = "cox", alpha = 0, nfolds = 5)
  model <- glmnet(X_train, Surv_train, family = "cox", alpha = 0, lambda = cv_ridge$lambda.min)
  
  p_val  <- predict(model, newx = X_val,  type = "link")
  p_test <- predict(model, newx = X_test, type = "link")
  
  c_val  <- concordance.index(p_val,  surv.time = time_val,  surv.event = event_val)$c.index
  c_test <- concordance.index(p_test, surv.time = time_test, surv.event = event_test)$c.index
  
  c(c_val, c_test)
}

k_vals <- seq(5, 51, by = 1)
for (k in k_vals) {
  top_genes <- head(rank1$Gene, k)
  cc <- fit_ridge_topk(top_genes, data1)
  all_topk_results <- rbind(all_topk_results, data.frame(
    Model = "Ridge", Dataset = "Ridge_51", k = k, Val_C = cc[1], Test_C = cc[2]
  ))
  cat(sprintf("k=%2d: Val=%.4f | Test=%.4f\n", k, cc[1], cc[2]))
}
cat("\n")


fit_lasso_topk <- function(genes, dat) {
  genes <- intersect(genes, colnames(dat))
  if (length(genes) == 0) return(c(NA, NA))
  
  X_train <- scale(as.matrix(dat[train_idx, genes, drop = FALSE]))
  X_val   <- scale(as.matrix(dat[val_idx,   genes, drop = FALSE]))
  X_test  <- scale(as.matrix(dat[test_idx,  genes, drop = FALSE]))
  
  cv_lasso <- cv.glmnet(X_train, Surv_train, family = "cox", alpha = 1, nfolds = 5)
  model <- glmnet(X_train, Surv_train, family = "cox", alpha = 1, lambda = cv_lasso$lambda.min)
  
  p_val  <- predict(model, newx = X_val,  type = "link")
  p_test <- predict(model, newx = X_test, type = "link")
  
  c_val  <- concordance.index(p_val,  surv.time = time_val,  surv.event = event_val)$c.index
  c_test <- concordance.index(p_test, surv.time = time_test, surv.event = event_test)$c.index
  
  c(c_val, c_test)
}

k_vals_lasso <- seq(2, min(60, nrow(rank2)), by = 1)
for (k in k_vals_lasso) {
  top_genes <- head(rank2$Gene, k)
  cc <- fit_lasso_topk(top_genes, data2)
  all_topk_results <- rbind(all_topk_results, data.frame(
    Model = "Lasso", Dataset = "Lasso_20", k = k, Val_C = cc[1], Test_C = cc[2]
  ))
  cat(sprintf("k=%2d: Val=%.4f | Test=%.4f\n", k, cc[1], cc[2]))
}
cat("\n")


fit_enet_topk <- function(genes, dat) {
  genes <- intersect(genes, colnames(dat))
  if (length(genes) == 0) return(c(NA, NA))
  
  X_train <- scale(as.matrix(dat[train_idx, genes, drop = FALSE]))
  X_val   <- scale(as.matrix(dat[val_idx,   genes, drop = FALSE]))
  X_test  <- scale(as.matrix(dat[test_idx,  genes, drop = FALSE]))
  
  cv_enet <- cv.glmnet(X_train, Surv_train, family = "cox", alpha = 0.5, nfolds = 5)
  model <- glmnet(X_train, Surv_train, family = "cox", alpha = 0.5, lambda = cv_enet$lambda.min)
  
  p_val  <- predict(model, newx = X_val,  type = "link")
  p_test <- predict(model, newx = X_test, type = "link")
  
  c_val  <- concordance.index(p_val,  surv.time = time_val,  surv.event = event_val)$c.index
  c_test <- concordance.index(p_test, surv.time = time_test, surv.event = event_test)$c.index
  
  c(c_val, c_test)
}

k_vals_enet <- seq(2, min(60, nrow(rank2)), by = 1)
for (k in k_vals_enet) {
  top_genes <- head(rank2$Gene, k)
  cc <- fit_enet_topk(top_genes, data2)
  all_topk_results <- rbind(all_topk_results, data.frame(
    Model = "ElasticNet", Dataset = "Lasso_20", k = k, Val_C = cc[1], Test_C = cc[2]
  ))
  cat(sprintf("k=%2d: Val=%.4f | Test=%.4f\n", k, cc[1], cc[2]))
}
cat("\n")

fit_stepwise_topk <- function(genes, dat) {
  genes <- intersect(genes, colnames(dat))
  if (length(genes) == 0) return(c(NA, NA))
  
  train_df <- data.frame(dat[train_idx, genes, drop = FALSE],
                         time = time_train, event = event_train)
  
  formula_str <- paste("Surv(time, event) ~", paste(genes, collapse = " + "))
  full_cox <- coxph(as.formula(formula_str), data = train_df)
  step_cox <- step(full_cox, direction = "forward", trace = 0)
  
  val_df  <- data.frame(dat[val_idx,  genes, drop = FALSE])
  test_df <- data.frame(dat[test_idx, genes, drop = FALSE])
  
  p_val  <- predict(step_cox, newdata = val_df,  type = "lp")
  p_test <- predict(step_cox, newdata = test_df, type = "lp")
  
  c_val  <- concordance.index(p_val,  surv.time = time_val,  surv.event = event_val)$c.index
  c_test <- concordance.index(p_test, surv.time = time_test, surv.event = event_test)$c.index
  
  c(c_val, c_test)
}

k_vals_step <- seq(5, 51, by = 2)
for (k in k_vals_step) {
  top_genes <- head(rank1$Gene, k)
  cc <- fit_stepwise_topk(top_genes, data1)
  all_topk_results <- rbind(all_topk_results, data.frame(
    Model = "StepwiseCox", Dataset = "Ridge_51", k = k, Val_C = cc[1], Test_C = cc[2]
  ))
  cat(sprintf("k=%2d: Val=%.4f | Test=%.4f\n", k, cc[1], cc[2]))
}
cat("\n")


fit_coxboost_topk <- function(genes, dat) {
  genes <- intersect(genes, colnames(dat))
  if (length(genes) == 0) return(c(NA, NA))
  
  X_train <- scale(as.matrix(dat[train_idx, genes, drop = FALSE]))
  X_val   <- scale(as.matrix(dat[val_idx,   genes, drop = FALSE]))
  X_test  <- scale(as.matrix(dat[test_idx,  genes, drop = FALSE]))
  
  tryCatch({
    cb_model <- CoxBoost(time = time_train, status = event_train,
                         x = X_train, stepno = 50, penalty = 100)
    p_val  <- as.vector(predict(cb_model, newdata = X_val,  at.step = 50))
    p_test <- as.vector(predict(cb_model, newdata = X_test, at.step = 50))
    
    keep_val  <- complete.cases(p_val,  time_val,  event_val)
    keep_test <- complete.cases(p_test, time_test, event_test)
    
    c_val  <- concordance.index(p_val[keep_val],
                                surv.time = time_val[keep_val],
                                surv.event = event_val[keep_val])$c.index
    c_test <- concordance.index(p_test[keep_test],
                                surv.time = time_test[keep_test],
                                surv.event = event_test[keep_test])$c.index
    
    c(c_val, c_test)
  }, error = function(e) {
    c(NA, NA)
  })
}

k_vals_cb <- seq(5, 41, by = 2)
for (k in k_vals_cb) {
  top_genes <- head(rank3$Gene, k)
  cc <- fit_coxboost_topk(top_genes, data3)
  all_topk_results <- rbind(all_topk_results, data.frame(
    Model = "CoxBoost", Dataset = "CoxBoost_41", k = k, Val_C = cc[1], Test_C = cc[2]
  ))
  cat(sprintf("k=%2d: Val=%.4f | Test=%.4f\n", k, cc[1], cc[2]))
}
cat("\n")


fit_rsf_topk <- function(genes, dat) {
  genes <- intersect(genes, colnames(dat))
  if (length(genes) == 0) return(c(NA, NA))
  
  train_df_rsf <- data.frame(dat[train_idx, genes, drop = FALSE],
                             time = time_train, event = event_train)
  
  tryCatch({
    rsf_model <- rfsrc(Surv(time, event) ~ .,
                       data = train_df_rsf,
                       ntree = 200,
                       nodesize = 5,
                       importance = TRUE)
    
    val_df_rsf  <- data.frame(dat[val_idx,  genes, drop = FALSE])
    test_df_rsf <- data.frame(dat[test_idx, genes, drop = FALSE])
    
    p_val  <- predict(rsf_model, newdata = val_df_rsf)$predicted
    p_test <- predict(rsf_model, newdata = test_df_rsf)$predicted
    
    c_val  <- concordance.index(p_val,  surv.time = time_val,  surv.event = event_val)$c.index
    c_test <- concordance.index(p_test, surv.time = time_test, surv.event = event_test)$c.index
    
    c(c_val, c_test)
  }, error = function(e) {
    c(NA, NA)
  })
}

k_vals_rsf <- seq(2, min(60, nrow(rank4)), by = 1)
for (k in k_vals_rsf) {
  top_genes <- head(rank4$Gene, k)
  cc <- fit_rsf_topk(top_genes, data4)
  all_topk_results <- rbind(all_topk_results, data.frame(
    Model = "RSF", Dataset = "RSF_20", k = k, Val_C = cc[1], Test_C = cc[2]
  ))
  cat(sprintf("k=%2d: Val=%.4f | Test=%.4f\n", k, cc[1], cc[2]))
}
cat("\n")


fit_plsr_topk <- function(genes, dat) {
  genes <- intersect(genes, colnames(dat))
  if (length(genes) == 0) return(c(NA, NA))
  
  X_train <- as.matrix(dat[train_idx, genes, drop = FALSE])
  X_val   <- as.matrix(dat[val_idx,   genes, drop = FALSE])
  X_test  <- as.matrix(dat[test_idx,  genes, drop = FALSE])
  
  tryCatch({
    plsr_model <- plsRcox(Xplan = X_train,
                          time = time_train,
                          event = event_train,
                          nt = min(3, length(genes)),
                          sparse = TRUE)
    
    p_val  <- predict(plsr_model, newdata = X_val)
    p_test <- predict(plsr_model, newdata = X_test)
    
    c_val  <- concordance.index(p_val,  surv.time = time_val,  surv.event = event_val)$c.index
    c_test <- concordance.index(p_test, surv.time = time_test, surv.event = event_test)$c.index
    
    c(c_val, c_test)
  }, error = function(e) {
    c(NA, NA)
  })
}

k_vals_pls <- seq(2, min(60, nrow(rank5)), by = 1)
for (k in k_vals_pls) {
  top_genes <- head(rank5$Gene, k)
  cc <- fit_plsr_topk(top_genes, data5)
  all_topk_results <- rbind(all_topk_results, data.frame(
    Model = "plsRcox", Dataset = "PLS_20", k = k, Val_C = cc[1], Test_C = cc[2]
  ))
  cat(sprintf("k=%2d: Val=%.4f | Test=%.4f\n", k, cc[1], cc[2]))
}
cat("\n")


fit_superpc_topk <- function(genes, dat) {
  genes <- intersect(genes, colnames(dat))
  if (length(genes) == 0) return(c(NA, NA))
  
  X_train <- as.matrix(dat[train_idx, genes, drop = FALSE])
  X_val   <- as.matrix(dat[val_idx,   genes, drop = FALSE])
  X_test  <- as.matrix(dat[test_idx,  genes, drop = FALSE])
  
  tryCatch({
    superpc_data <- list(
      x = t(X_train),
      y = time_train,
      censoring.status = event_train,
      featurenames = colnames(X_train)
    )
    
    sp_train <- superpc.train(superpc_data, type = "survival")
    
    sp_val_data <- list(
      x = t(X_val),
      y = time_val,
      censoring.status = event_val,
      featurenames = colnames(X_val)
    )
    
    sp_test_data <- list(
      x = t(X_test),
      y = time_test,
      censoring.status = event_test,
      featurenames = colnames(X_test)
    )
    
    sp_pred_val  <- superpc.predict(sp_train, superpc_data, sp_val_data,
                                    threshold = 1.5, n.components = 1)
    sp_pred_test <- superpc.predict(sp_train, superpc_data, sp_test_data,
                                    threshold = 1.5, n.components = 1)
    
    c_val  <- concordance.index(sp_pred_val$v.pred,
                                surv.time = time_val,
                                surv.event = event_val)$c.index
    c_test <- concordance.index(sp_pred_test$v.pred,
                                surv.time = time_test,
                                surv.event = event_test)$c.index
    
    c(c_val, c_test)
  }, error = function(e) {
    c(NA, NA)
  })
}

k_vals_sp <- seq(2, min(60, nrow(rank5)), by = 1)
for (k in k_vals_sp) {
  top_genes <- head(rank5$Gene, k)
  cc <- fit_superpc_topk(top_genes, data5)
  all_topk_results <- rbind(all_topk_results, data.frame(
    Model = "SuperPC", Dataset = "PLS_20", k = k, Val_C = cc[1], Test_C = cc[2]
  ))
  cat(sprintf("k=%2d: Val=%.4f | Test=%.4f\n", k, cc[1], cc[2]))
}
cat("\n")


fit_gbm_topk <- function(genes, dat) {
  genes <- intersect(genes, colnames(dat))
  if (length(genes) == 0) return(c(NA, NA))
  
  train_df_gbm <- data.frame(dat[train_idx, genes, drop = FALSE])
  train_df_gbm$time  <- time_train
  train_df_gbm$event <- event_train
  
  tryCatch({
    gbm_model <- gbm(
      formula = Surv(time, event) ~ .,
      data = train_df_gbm,
      distribution = "coxph",
      n.trees = 300,
      interaction.depth = 3,
      shrinkage = 0.01,
      n.minobsinnode = 10,
      verbose = FALSE
    )
    
    best_trees <- gbm.perf(gbm_model, method = "OOB", plot.it = FALSE)
    
    val_df_gbm  <- data.frame(dat[val_idx,  genes, drop = FALSE])
    test_df_gbm <- data.frame(dat[test_idx, genes, drop = FALSE])
    
    p_val  <- predict(gbm_model, newdata = val_df_gbm,
                      n.trees = best_trees, type = "response")
    p_test <- predict(gbm_model, newdata = test_df_gbm,
                      n.trees = best_trees, type = "response")
    
    c_val  <- concordance.index(p_val,  surv.time = time_val,  surv.event = event_val)$c.index
    c_test <- concordance.index(p_test, surv.time = time_test, surv.event = event_test)$c.index
    
    c(c_val, c_test)
  }, error = function(e) {
    c(NA, NA)
  })
}

k_vals_gbm <- seq(2, min(60, nrow(rank4)), by = 1)
for (k in k_vals_gbm) {
  top_genes <- head(rank4$Gene, k)
  cc <- fit_gbm_topk(top_genes, data4)
  all_topk_results <- rbind(all_topk_results, data.frame(
    Model = "GBM", Dataset = "RSF_20", k = k, Val_C = cc[1], Test_C = cc[2]
  ))
  cat(sprintf("k=%2d: Val=%.4f | Test=%.4f\n", k, cc[1], cc[2]))
}
cat("\n")


print(all_topk_results)


best_per_model <- all_topk_results %>%
  group_by(Model) %>%
  slice_max(Test_C, n = 1) %>%
  ungroup() %>%
  arrange(desc(Test_C))

print(best_per_model)
write.csv(all_topk_results, "topk_optimization_all_models.csv", row.names = FALSE)
write.csv(best_per_model, "/path/to/topk_optimization_best_k_per_model.csv", row.names = FALSE)

