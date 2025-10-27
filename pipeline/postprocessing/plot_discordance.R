library(ape)
library(dplyr)
library(stringr)

source("pipeline/plot_config.R")
source("src/discordance.R")

regressand <- "discordance_aa_rbd_vs_ntd"
model_path_aa <- file.path(brmsdir, "models", paste0(regressand, ".Rds"))
model_aa <- if (file.exists(model_path_aa)) readRDS(model_path_aa) else NULL
out_aa <- discordance_analysis(
  rbd_tree_aa,
  ntd_tree_aa,
  virusinfo,
  regressand = regressand,
  model = model_aa,
  signif_only = TRUE
)

regressand <- "discordance_3di_rbd_vs_ntd"
model_path_3di <- file.path(brmsdir, "models", paste0(regressand, ".Rds"))
model_3di <- if (file.exists(model_path_3di)) readRDS(model_path_3di) else NULL
out_3di <- discordance_analysis(
  rbd_tree_3di,
  ntd_tree_3di,
  virusinfo,
  model = model_3di,
  regressand = regressand,
  signif_only = TRUE
)

df_aa <- out_aa$df %>% mutate(dataset = "aa")
df_3di <- out_3di$df %>% mutate(dataset = "3di")
combined_df <- bind_rows(df_aa, df_3di)
combined_discordance_plot(df_aa, file.path(imgdir, "combined_discordance_plot_aa.pdf"))
combined_discordance_plot(df_3di, file.path(imgdir, "combined_discordance_plot_3di.pdf"))

cat("Plotting combined estimates and correlations...\n")
combined_estimate_plot(
  out_aa$model, out_3di$model,
  filename = file.path(imgdir, "combined_estimate_plot.pdf"),
  with_intercept = FALSE
)

combined_corr_plot(
  out_aa$cormat, out_3di$cormat,
  file.path(imgdir, "combined_corr_plot.pdf")
)
