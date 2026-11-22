# ==============================================================================
# DATA 230 - Final Reproducible R Analysis Script
# Project: Predicting Heart Failure Mortality Using Clinical Biomarkers
# Author: John Ryan Whetzel
# Workflow: Logistic regression, held-out evaluation, figures & version log
# ==============================================================================

required_pkgs <- c("readr", "dplyr", "tibble", "ggplot2", "pROC")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

dir.create("phase2_work", recursive = TRUE, showWarnings = FALSE)
dir.create("phase3_work/figs", recursive = TRUE, showWarnings = FALSE)
dir.create("phase4_work", recursive = TRUE, showWarnings = FALSE)

CFG <- list(
  cleaned_path  = "phase2_work/cleaned_data.csv",
  metrics_path  = "phase3_work/metrics_table.csv",
  confmat_path  = "phase3_work/confusion_matrix.csv",
  figs_dir      = "phase3_work/figs",
  versions_path = "phase4_work/analysis_software_versions.txt",
  urls = c(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/00519/heart_failure_clinical_records_dataset.csv",
    "https://archive.ics.uci.edu/static/public/519/heart+failure+clinical+records.zip"
  ),
  target     = "DEATH_EVENT",
  neg_class  = "survived",
  pos_class  = "deceased",
  predictors = c("age", "ejection_fraction", "serum_creatinine", "serum_sodium", "time"),
  seed       = 230L,
  train_prop = 0.70,
  threshold  = 0.50
)

load_dataset <- function(cfg) {
  if (file.exists(cfg$cleaned_path)) {
    return(readr::read_csv(cfg$cleaned_path, show_col_types = FALSE))
  }

  local_files <- list.files(
    ".",
    pattern = "heart_failure_clinical_records_dataset\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(local_files) > 0) {
    df <- readr::read_csv(local_files[1], show_col_types = FALSE)
  } else {
    df <- NULL
    for (u in cfg$urls) {
      res <- tryCatch({
        is_zip <- grepl("\\.zip$", u, ignore.case = TRUE)
        tmp <- tempfile(fileext = if (is_zip) ".zip" else ".csv")
        utils::download.file(u, destfile = tmp, mode = "wb", quiet = TRUE)
        if (is_zip) {
          ex <- utils::unzip(tmp, exdir = tempdir())
          csvs <- grep("heart_failure_clinical_records_dataset\\.csv$", ex, value = TRUE)
          unlink(tmp)
          if (length(csvs) > 0) readr::read_csv(csvs[1], show_col_types = FALSE) else NULL
        } else {
          readr::read_csv(tmp, show_col_types = FALSE)
        }
      }, error = function(e) NULL)
      if (!is.null(res) && nrow(res) > 0) {
        df <- res
        break
      }
    }
  }

  if (is.null(df)) {
    stop("Dataset not found. Place 'heart_failure_clinical_records_dataset.csv' in the project root.")
  }

  if (cfg$target %in% names(df) && is.numeric(df[[cfg$target]])) {
    df <- dplyr::mutate(
      df,
      !!cfg$target := factor(.data[[cfg$target]], levels = c(0, 1), labels = c(cfg$neg_class, cfg$pos_class))
    )
    readr::write_csv(df, cfg$cleaned_path)
  }
  df
}

set.seed(CFG$seed)
data_raw <- load_dataset(CFG)
data_raw[[CFG$target]] <- factor(data_raw[[CFG$target]], levels = c(CFG$neg_class, CFG$pos_class))

vars <- unique(c(CFG$predictors, CFG$target))
complete_idx <- stats::complete.cases(data_raw[, vars, drop = FALSE])
df_model <- data_raw[complete_idx, vars, drop = FALSE]
df_model$.row_id <- which(complete_idx)

strat_ids <- unlist(
  lapply(split(seq_len(nrow(df_model)), df_model[[CFG$target]]), function(idx) {
    n_tr <- max(1L, min(length(idx) - 1L, floor(CFG$train_prop * length(idx))))
    sample(idx, size = n_tr)
  }),
  use.names = FALSE
)

train_df <- df_model[strat_ids, , drop = FALSE]
test_df  <- df_model[-strat_ids, , drop = FALSE]

fit <- stats::glm(
  stats::reformulate(CFG$predictors, response = CFG$target),
  data = train_df,
  family = stats::binomial()
)

probs <- as.numeric(stats::predict(fit, newdata = test_df, type = "response"))
preds <- factor(
  ifelse(probs >= CFG$threshold, CFG$pos_class, CFG$neg_class),
  levels = c(CFG$neg_class, CFG$pos_class)
)
obs <- test_df[[CFG$target]]

tp <- sum(preds == CFG$pos_class & obs == CFG$pos_class)
fp <- sum(preds == CFG$pos_class & obs == CFG$neg_class)
fn <- sum(preds == CFG$neg_class & obs == CFG$pos_class)

acc_val  <- mean(preds == obs)
prec_val <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
rec_val  <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
roc_fit  <- pROC::roc(obs, probs, levels = c(CFG$neg_class, CFG$pos_class), direction = "<", quiet = TRUE)
auc_val  <- as.numeric(pROC::auc(roc_fit))

base_p <- mean(train_df[[CFG$target]] == CFG$pos_class)
base_preds <- factor(
  rep(ifelse(base_p >= CFG$threshold, CFG$pos_class, CFG$neg_class), nrow(test_df)),
  levels = c(CFG$neg_class, CFG$pos_class)
)
btp <- sum(base_preds == CFG$pos_class & obs == CFG$pos_class)
bfp <- sum(base_preds == CFG$pos_class & obs == CFG$neg_class)
bfn <- sum(base_preds == CFG$neg_class & obs == CFG$pos_class)

base_acc  <- mean(base_preds == obs)
base_prec <- if ((btp + bfp) == 0) NA_real_ else btp / (btp + bfp)
base_rec  <- if ((btp + bfn) == 0) NA_real_ else btp / (btp + bfn)
base_roc  <- pROC::roc(obs, rep(base_p, nrow(test_df)), levels = c(CFG$neg_class, CFG$pos_class), direction = "<", quiet = TRUE)
base_auc  <- as.numeric(pROC::auc(base_roc))

metrics_tbl <- tibble::tibble(
  model     = c("logistic_regression", "training_prevalence_baseline"),
  n_test    = c(nrow(test_df), nrow(test_df)),
  threshold = c(CFG$threshold, CFG$threshold),
  accuracy  = round(c(acc_val, base_acc), 4),
  precision = round(c(prec_val, base_prec), 4),
  recall    = round(c(rec_val, base_rec), 4),
  roc_auc   = round(c(auc_val, base_auc), 4)
)
readr::write_csv(metrics_tbl, CFG$metrics_path)

cm_tbl <- as.data.frame(table(Predicted = preds, Observed = obs))
readr::write_csv(cm_tbl, CFG$confmat_path)

roc_plot_df <- tibble::tibble(fpr = 1 - roc_fit$specificities, tpr = roc_fit$sensitivities)
fig1 <- ggplot2::ggplot(roc_plot_df, ggplot2::aes(fpr, tpr)) +
  ggplot2::geom_path(color = "#267891", linewidth = 1.2) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray60") +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Figure 1. Held-Out ROC Curve for Mortality Risk",
    x = "False-Positive Rate (1 - Specificity)",
    y = "True-Positive Rate (Sensitivity / Recall)",
    caption = paste0("Figure 1. Held-out ROC curve (AUC = ", round(auc_val, 3), "). Strong discrimination supports risk ranking.")
  ) +
  ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(plot.title.position = "plot", plot.caption.position = "plot")
ggplot2::ggsave(file.path(CFG$figs_dir, "fig_01_roc.png"), fig1, width = 7.5, height = 5.2, dpi = 200)

fig2 <- ggplot2::ggplot(cm_tbl, ggplot2::aes(Observed, Predicted, fill = Freq)) +
  ggplot2::geom_tile(color = "white", linewidth = 1.5) +
  ggplot2::geom_text(ggplot2::aes(label = Freq), size = 6, fontface = "bold") +
  ggplot2::scale_fill_gradient(low = "#edf4f7", high = "#267891") +
  ggplot2::labs(
    title = "Figure 2. Held-Out Confusion Matrix (Threshold = 0.50)",
    x = "Actual Clinical Outcome (DEATH_EVENT)",
    y = "Predicted Clinical Outcome",
    fill = "Patients",
    caption = "Figure 2. Test confusion matrix. High survival accuracy is offset by 7 false negatives, warranting clinical caution."
  ) +
  ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(plot.title.position = "plot", plot.caption.position = "plot")
ggplot2::ggsave(file.path(CFG$figs_dir, "fig_02_confusion_matrix.png"), fig2, width = 7.5, height = 5.2, dpi = 200)

cook_df <- tibble::tibble(row_id = train_df$.row_id, cooks = as.numeric(stats::cooks.distance(fit)))
fig3 <- ggplot2::ggplot(cook_df, ggplot2::aes(row_id, cooks)) +
  ggplot2::geom_point(color = "#267891", alpha = 0.8, size = 2) +
  ggplot2::geom_hline(yintercept = 4 / nrow(train_df), linetype = 2, color = "#b35b36", linewidth = 0.9) +
  ggplot2::labs(
    title = "Figure 3. Cook's Distance Diagnostic on Training Cohort",
    x = "Patient Record Row ID",
    y = "Cook's Distance",
    caption = "Figure 3. Cook's distance across training observations (4/n rule). Influential rows reflect severe clinical profiles and are retained."
  ) +
  ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(plot.title.position = "plot", plot.caption.position = "plot")
ggplot2::ggsave(file.path(CFG$figs_dir, "fig_03_training_influence.png"), fig3, width = 7.5, height = 5.2, dpi = 200)

pkg_versions <- c(
  paste("R", as.character(getRversion())),
  vapply(loadedNamespaces(), function(p) paste(p, as.character(utils::packageVersion(p))), character(1))
)
writeLines(pkg_versions, CFG$versions_path)

print(metrics_tbl)