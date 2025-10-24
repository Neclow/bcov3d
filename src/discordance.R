library(brms)
library(dplyr)
library(TreeDist)

source("src/plot.R")
source("src/treestats.R")

merge_discordance_with_tree_stats <- function(tree1, tree2, df) {
  pathnode1 <- as.data.frame(pathnode(tree1))
  rownames(pathnode1) <- tree1$tip.label
  pathnode2 <- as.data.frame(pathnode(tree2))
  rownames(pathnode2) <- tree2$tip.label

  pathnodes <- merge(pathnode1, pathnode2, by = "row.names", suffixes = c("_1", "_2"))
  merged_df <- merge(
    x = df,
    y = pathnodes,
    by.x = "Label",
    by.y = "Row.names"
  )
  merged_df$roottotip_avg <- (merged_df$roottotippath_1 + merged_df$roottotippath_2) / 2
  merged_df %>%
    mutate(across(
      c(Contribution, roottotip_avg),
      ~ as.numeric(scale(.x))
    )) %>%
    mutate(is_sarb = Subgenus == "Sarb.")
}

get_discordance_per_taxon <- function(tree1, tree2) {
  taxa <- intersect(tree1$tip.label, tree2$tip.label)
  taxa_discordances <- vector()
  for (i in seq_along(taxa)) {
    mci <- MutualClusteringInfo(
      drop.tip(tree1, taxa[i]),
      drop.tip(tree2, taxa[i])
    )
    taxa_discordances[i] <- mci
  }
  names(taxa_discordances) <- taxa
  taxa_discordances
}

group_discordance_by_receptor <- function(tree1, tree2, virusinfo) {
  contrib <- get_discordance_per_taxon(tree1, tree2)
  data <- list()
  for (i in seq_along(tree1$tip.label)) {
    label <- tree1$tip.label[[i]]
    parts <- strsplit(label, "/")[[1]]
    virus_name <- gsub(" ", "", paste(head(parts, -2), collapse = "/"))
    receptor <- virusinfo %>%
      filter(Virus == virus_name) %>%
      select(Receptor) %>%
      slice(1) %>%
      pull()
    subgenus <- virusinfo %>%
      filter(Virus == virus_name) %>%
      select(Subgenus) %>%
      slice(1) %>%
      pull()
    data[[i]] <- c(label, receptor, subgenus, contrib[[label]])
  }
  df <- as.data.frame(do.call(rbind, data))
  colnames(df) <- c("Label", "Receptor", "Subgenus", "Contribution")
  df$Contribution <- as.numeric(df$Contribution)
  df %>%
    mutate(is_ace2 = Receptor == "ACE2")
}

brms_save <- function(model, regressand) {
  models_dir <- file.path(brmsdir, "models")
  summary_dir <- file.path(brmsdir, "summary")
  preds_dir <- file.path(brmsdir, "preds")
  dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(preds_dir, showWarnings = FALSE, recursive = TRUE)
  # Save model
  saveRDS(
    model,
    paste0(models_dir, "/", regressand, ".Rds", sep = "")
  )
  # Save summary and predictions
  write.csv(
    summary(model)$fixed,
    paste0(summary_dir, "/", regressand, ".csv", sep = "")
  )
  write.csv(
    predict(model),
    paste0(preds_dir, "/", regressand, ".csv", sep = "")
  )
}

discordance_analysis <- function(
  tree1, tree2, virusinfo, regressand, model = NULL, signif_only = FALSE
) {
  df <- group_discordance_by_receptor(tree1, tree2, virusinfo)
  cat("Plotting discordance...\n")
  discordance_plot(df, file.path(imgdir, paste0(regressand, ".pdf")))
  merged_df <- merge_discordance_with_tree_stats(tree1, tree2, df)

  cor_df <- cor(
    merged_df[, c("roottotippath_1", "roottotippath_2", "Contribution")]
  )

  corr_plot(cor_df, regressand)
  
  if (is.null(model)) {
    cat("Modelling with BRMS...\n")
    formula <- bf(
      is_ace2 ~ Contribution + roottotippath_1 + roottotippath_2,
      family = bernoulli(link = "logit")
    )
    model <- brm(formula,
      data = merged_df,
      chains = 4,
      iter = 4000,
      warmup = 2000,
      control = list(
        adapt_delta = 0.99,
        max_treedepth = 15
      ),
      cores = 4
    )
  }
  cat("Saving brms outputs...\n")
  brms_save(model, regressand)
  cat("Plotting brms outputs...\n")
  brms_plot(model, merged_df, regressand, signif_only = signif_only)

  list(df = merged_df, cormat = cor_df, model = model)
}
