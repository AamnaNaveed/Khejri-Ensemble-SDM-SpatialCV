# ============================================
# ENSEMBLE SDM WITH SPATIAL BLOCK CV
# For: Prosopis cineraria in Pakistan
# Author: Aamna Naveed
# ============================================

library(dismo)
library(raster)
library(ggplot2)
library(sf)
library(randomForest)

setwd("C:/Users/Laptop Land/Desktop/Khejri-Ensemble-SDM-Advanced")

# --- CONFIGURATION ---
selected <- read.csv("selected_variables.csv")$Variable

# --- LOAD DATA ---
presence <- read.csv("khejri_presence_with_climate.csv")
background <- read.csv("khejri_background_with_climate.csv")

pres <- presence[
  presence$decimalLongitude >= 60 & presence$decimalLongitude <= 78 &
  presence$decimalLatitude >= 23 & presence$decimalLatitude <= 38,
]

pres_coords <- pres[, c("decimalLongitude", "decimalLatitude")]
pres_vars <- pres[, selected]
pres_vars <- na.omit(pres_vars)
pres_coords <- pres_coords[rownames(pres_vars), ]

bg_coords <- background[, c("decimalLongitude", "decimalLatitude")]
bg_vars <- background[, selected]
bg_vars <- na.omit(bg_vars)
bg_coords <- bg_coords[rownames(bg_vars), ]

# --- SPATIAL BLOCKS (4 folds by longitude) ---
assign_block <- function(lon) {
  ifelse(lon < 64.5, 1, ifelse(lon < 69, 2, ifelse(lon < 73.5, 3, 4)))
}
pres_blocks <- assign_block(pres_coords[, "decimalLongitude"])
bg_blocks <- assign_block(bg_coords[, "decimalLongitude"])

# --- CROSS-VALIDATION FUNCTION ---
run_fold <- function(fold_num, pres_vars, pres_coords, bg_vars, bg_coords, 
                     pres_blocks, bg_blocks, selected) {
  
  train_pres_idx <- which(pres_blocks != fold_num)
  test_pres_idx <- which(pres_blocks == fold_num)
  train_bg_idx <- which(bg_blocks != fold_num)
  test_bg_idx <- which(bg_blocks == fold_num)
  
  if (length(test_pres_idx) < 3 || length(test_bg_idx) < 10) return(NULL)
  
  pres_train <- pres_vars[train_pres_idx, ]
  pres_test <- pres_vars[test_pres_idx, ]
  bg_train <- bg_vars[train_bg_idx, ]
  bg_test <- bg_vars[test_bg_idx, ]
  
  # BIOCLIM
  bc <- bioclim(pres_train)
  bc_pred_pres <- predict(bc, pres_test)
  bc_pred_bg <- predict(bc, bg_test)
  bc_min <- min(c(bc_pred_pres, bc_pred_bg), na.rm = TRUE)
  bc_max <- max(c(bc_pred_pres, bc_pred_bg), na.rm = TRUE)
  bc_pred_pres_std <- (bc_pred_pres - bc_min) / (bc_max - bc_min)
  bc_pred_bg_std <- (bc_pred_bg - bc_min) / (bc_max - bc_min)
  
  # GLM
  train_glm <- rbind(pres_train, bg_train)
  train_glm$pa <- c(rep(1, nrow(pres_train)), rep(0, nrow(bg_train)))
  glm_mod <- glm(pa ~ ., data = train_glm, family = binomial)
  glm_pred_pres <- predict(glm_mod, pres_test, type = "response")
  glm_pred_bg <- predict(glm_mod, bg_test, type = "response")
  
  # Random Forest
  train_rf <- train_glm
  train_rf$pa <- as.factor(train_rf$pa)
  rf_mod <- randomForest(pa ~ ., data = train_rf, ntree = 500)
  rf_pred_pres <- predict(rf_mod, pres_test, type = "prob")[, "1"]
  rf_pred_bg <- predict(rf_mod, bg_test, type = "prob")[, "1"]
  
  # Evaluation
  e_bc <- evaluate(p = bc_pred_pres_std, a = bc_pred_bg_std)
  e_glm <- evaluate(p = glm_pred_pres, a = glm_pred_bg)
  e_rf <- evaluate(p = rf_pred_pres, a = rf_pred_bg)
  
  calc_tss <- function(eval_obj) {
    cm <- eval_obj@confusion
    tpr <- cm[, "tp"] / (cm[, "tp"] + cm[, "fn"])
    tnr <- cm[, "tn"] / (cm[, "tn"] + cm[, "fp"])
    return(max(tpr + tnr - 1, na.rm = TRUE))
  }
  
  # Equal-weight ensemble
  ens_pred_pres <- (bc_pred_pres_std + glm_pred_pres + rf_pred_pres) / 3
  ens_pred_bg <- (bc_pred_bg_std + glm_pred_bg + rf_pred_bg) / 3
  e_ens <- evaluate(p = ens_pred_pres, a = ens_pred_bg)
  
  return(list(
    auc_bc = e_bc@auc, auc_glm = e_glm@auc, auc_rf = e_rf@auc, auc_ens = e_ens@auc,
    tss_bc = calc_tss(e_bc), tss_glm = calc_tss(e_glm), 
    tss_rf = calc_tss(e_rf), tss_ens = calc_tss(e_ens)
  ))
}

# --- RUN CV ---
results_list <- list()
for (fold in 1:4) {
  res <- run_fold(fold, pres_vars, pres_coords, bg_vars, bg_coords,
                  pres_blocks, bg_blocks, selected)
  if (!is.null(res)) results_list[[length(results_list) + 1]] <- res
}

# --- COMPILE RESULTS ---
n_folds <- length(results_list)
results_df <- data.frame(
  Algorithm = c("BIOCLIM", "GLM", "Random Forest", "Ensemble"),
  AUC_Mean = round(c(mean(sapply(results_list, function(x) x$auc_bc)),
                     mean(sapply(results_list, function(x) x$auc_glm)),
                     mean(sapply(results_list, function(x) x$auc_rf)),
                     mean(sapply(results_list, function(x) x$auc_ens))), 3),
  AUC_SD = round(c(sd(sapply(results_list, function(x) x$auc_bc)),
                   sd(sapply(results_list, function(x) x$auc_glm)),
                   sd(sapply(results_list, function(x) x$auc_rf)),
                   sd(sapply(results_list, function(x) x$auc_ens))), 3)
)

print(results_df)
write.csv(results_df, "spatial_cv_results.csv", row.names = FALSE)

# --- BOXPLOT ---
plot_df <- data.frame(
  Algorithm = rep(c("BIOCLIM", "GLM", "Random Forest", "Ensemble"), each = n_folds),
  AUC = c(sapply(results_list, function(x) x$auc_bc),
          sapply(results_list, function(x) x$auc_glm),
          sapply(results_list, function(x) x$auc_rf),
          sapply(results_list, function(x) x$auc_ens))
)

p <- ggplot(plot_df, aes(x = reorder(Algorithm, AUC, median), y = AUC)) +
  geom_boxplot(aes(fill = Algorithm), alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(title = "Spatial Block Cross-Validation: Algorithm Performance (4 folds)",
       x = "Algorithm", y = "AUC") +
  theme_minimal() + theme(legend.position = "none")

ggsave("fig_01_spatial_cv_boxplot.png", p, width = 10, height = 7, dpi = 300)
cat("Done! Check spatial_cv_results.csv and fig_01_spatial_cv_boxplot.png\n")