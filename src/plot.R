library(brms)
library(dplyr)
library(envalysis)
library(ggplot2)
library(ggsci)
library(patchwork)
library(purrr)
library(reshape2)
library(tidybayes)

source("src/config.R")

extract_fixed <- function(model, with_intercept = TRUE) {
  fixed <- as.data.frame(summary(model)$fixed) %>%
    rename(l95 = "l-95% CI", u95 = "u-95% CI")
  if (!with_intercept) {
    fixed <- fixed %>% slice(-1)
  }
  fixed$term <- factor(rownames(fixed), levels = rev(rownames(fixed)))
  fixed
}

brms_plot <- function(model, data, regressand, signif_only = FALSE) {
  # BRMS convergence plot
  cat("Plot convergence...\n")
  pdf(file.path(imgdir, paste("brms_", regressand, ".pdf", sep = "")))
  plot(model, ask = FALSE)
  dev.off()

  # BRMS posterior predictive check
  cat("Plot posterior predictive check...\n")
  pdf(file.path(imgdir, paste("brms_ppc_", regressand, ".pdf", sep = "")),
    width = 2.63,
    height = 2.63
  )
  print(pp_check(model) + theme_publish())
  dev.off()

  # Predictions
  pdf(file.path(imgdir, paste("brms_epred_", regressand, ".pdf", sep = "")))
  cat("Plot predicted probabilities...\n")
  print(
    data %>%
      add_epred_draws(model) %>%
      ggplot(aes(x = .epred, fill = factor(is_ace2))) +
      geom_density(alpha = 0.7) +
      labs(
        x = "Predicted Probability", fill = "Observed is_ace2"
      )
  )
  dev.off()

  # # HDI
  # pdf(file.path(imgdir, paste("brms_hdi_", regressand, ".pdf", sep = "")))
  # hdi_range <- bayestestR::hdi(model, ci = c(0.5, 0.75, 0.89, 0.95))
  # print(plot(hdi_range, show_intercept = TRUE) + theme_publish() + scale_fill_npg())
  # dev.off()

  # Plot regressands
  cat("Plot regression coefficients...\n")
  pdf(file.path(imgdir, paste("brms_summary_", regressand, ".pdf", sep = "")),
    width = 4.63,
    height = 4.63
  )
  fixed <- extract_fixed(model)

  # Plot coefficients with error bars
  print(
    ggplot(fixed, aes(y = term, x = Estimate, color = term)) +
      geom_point(size = 3) +
      geom_errorbar(aes(xmin = l95, xmax = u95, color = term), width = 0.2) +
      labs(
        y = "Coefficient",
        x = "Estimate"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_publish() +
      scale_color_ucscgb()
  )
  dev.off()

  # Conditional effects
  model$data$is_ace2 <- as.numeric(model$data$is_ace2)
  if (signif_only) {
    effects <- c(as.data.frame(summary(model)$fixed) |>
      rename(l95 = "l-95% CI", u95 = "u-95% CI") |>
      filter((l95 < 0 & u95 < 0) | (l95 > 0 & u95 > 0)) |>
      rownames())
    # Remove the intercept if present
    effects <- effects[effects != "Intercept"]
  } else {
    effects <- NULL
  }

  #   map2(rev(pal_ucscgb("default")(2)), ~ .x +
  #     geom_smooth(color = .y, se = TRUE) +
  #     scale_color_ucscgb(alpha = 0.5))
  # if (is.null(effects) || length(effects) > 1) {
  #   p <- p %>%
  #     patchwork::wrap_plots(ncol = 1, axis_titles = "collect")
  # }
  cat("Plot conditional effects...\n")
  pdf(
    file.path(imgdir, paste("brms_effects_", regressand, ".pdf", sep = "")),
    width = 2.63, height = 2.63
  )
  p <- conditional_effects(model, prob = 0.95, effects = effects, plot = FALSE)
  print(length(effects))
  # Handle single vs multiple plots differently
  if (!is.null(effects) && length(effects) == 1) {
    # Single effect: p is already a single plot, just add styling
    p <- plot(
      p,
      points = TRUE,
      offset = TRUE,
      ask = FALSE,
      theme = theme_publish(),
      line_args = list(color = rev(pal_ucscgb("default", alpha = 0.25)(1))[1])
    )[[1]] +
      labs(y = "Uses ACE2 [Yes (1)/No (0)]") +
      labs(x = case_when(
        grepl("roottotippath_2", p$labels$x) ~ "RTT length (NTD)",
        grepl("roottotippath_1", p$labels$x) ~ "RTT length (RBD)",
        grepl("Contribution", p$labels$x) ~ "TDI"
      ))
  } else {
    # Multiple effects: apply colors and combine
    p <- plot(p,
      points = TRUE,
      offset = TRUE,
      ask = FALSE,
      theme = theme_publish()
    ) %>%
      map2(rev(pal_ucscgb("default")(length(p))), ~ .x +
        scale_color_ucscgb(alpha = 0.25) +
        labs(
          y = "Uses ACE2 [Yes (1)/No (0)]",
          x = case_when(
            grepl("roottotippath_2", .x$labels$x) ~ "RTT length (NTD)",
            grepl("roottotippath_1", .x$labels$x) ~ "RTT length (RBD)",
            grepl("Contribution", .x$labels$x) ~ "TDI",
            TRUE ~ .x$labels$x
          )
        )) %>%
      patchwork::wrap_plots(ncol = 1, axis_titles = "keep")
  }
  dev.off()
}

discordance_plot <- function(df, filename) {
  plot_data <- df %>%
    group_by(Receptor) %>%
    summarize(
      mean_contribution = mean(Contribution, na.rm = TRUE),
      sd_contribution = sd(Contribution, na.rm = TRUE),
      n = n(),
      # Calculate the standard error of the mean
      se = sd_contribution / sqrt(n),
      # Calculate the lower and upper bounds of the 95% CI
      lower_ci = mean_contribution - 1.96 * se,
      upper_ci = mean_contribution + 1.96 * se
    ) %>%
    ungroup() %>%
    # Sort the data by the mean_contribution in ascending order
    arrange(mean_contribution)

  # Factor for x-axis order
  plot_data$Receptor_x <- factor(plot_data$Receptor, levels = plot_data$Receptor)

  # Factor for color mapping (alphabetical)
  plot_data$Receptor_color <- factor(
    plot_data$Receptor,
    levels = sort(unique(plot_data$Receptor))
  )

  contrib_plot <- ggplot(
    plot_data,
    aes(x = Receptor_x, y = mean_contribution, color = Receptor_color)
  ) +
    # Add the points for the mean contribution
    geom_point(size = 4) +
    # Add the error bars for the 95% confidence interval
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    # Use ggsci's Nature Publishing Group color palette
    scale_color_npg() +
    # Customize the plot appearance
    labs(
      x = "Receptor",
      y = "Mean Contribution"
    ) +
    theme_publish() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave(
    filename = filename,
    plot = contrib_plot,
    device = "pdf",
    width = 4.63,
    height = 4.63
  )
}

combined_discordance_plot <- function(combined_df, filename) {
  plot_data <- combined_df %>%
    group_by(Receptor, dataset) %>% # Add dataset to grouping
    summarize(
      mean_contribution = mean(Contribution, na.rm = TRUE),
      sd_contribution = sd(Contribution, na.rm = TRUE),
      n = n(),
      se = sd_contribution / sqrt(n),
      lower_ci = mean_contribution - 1.96 * se,
      upper_ci = mean_contribution + 1.96 * se,
      .groups = "drop" # Ungroup automatically
    ) %>%
    # Sort by overall mean contribution across datasets for consistent x-axis ordering
    arrange(Receptor, dataset)

  # Create ordering based on overall mean contribution per receptor
  receptor_order <- combined_df %>%
    group_by(Receptor) %>%
    summarize(overall_mean = mean(Contribution, na.rm = TRUE), .groups = "drop") %>%
    arrange(overall_mean) %>%
    pull(Receptor)

  # Factor for x-axis order (based on overall mean)
  plot_data$Receptor_x <- factor(plot_data$Receptor, levels = receptor_order)
  # Factor for color mapping (alphabetical)
  plot_data$Receptor_color <- factor(
    plot_data$Receptor,
    levels = sort(unique(plot_data$Receptor))
  )

  contrib_plot <- ggplot(plot_data, aes(
    x = Receptor_x, y = mean_contribution,
    color = Receptor_color, alpha = dataset
  )) +
    geom_point(size = 4, position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
      width = 0.2, position = position_dodge(width = 0.3)
    ) +
    scale_color_npg() +
    scale_alpha_manual(values = c("aa" = 1.0, "3di" = 0.5)) +
    labs(
      x = "Receptor",
      y = "Mean Contribution"
    ) +
    theme_publish() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(
      alpha = guide_legend(title = "Dataset"),
      color = "none"
    )
  ggsave(
    filename = filename,
    plot = contrib_plot,
    device = "pdf",
    width = 4.63,
    height = 4.63
  )
}

corr_plot <- function(cormat, regressand) {
  cor_matrix <- matrix(cormat,
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("RTT length (RBD)", "RTT length (NTD)", "TDI"),
      c("RTT length (RBD)", "RTT length (NTD)", "TDI")
    )
  )

  # Get lower triangle only
  cor_matrix[upper.tri(cor_matrix)] <- NA

  # Convert to long format for ggplot2
  cor_melt <- melt(cor_matrix)
  colnames(cor_melt) <- c("Var1", "Var2", "value")

  # Remove NA values (upper triangle)
  cor_melt <- cor_melt[!is.na(cor_melt$value), ]

  pdf(file.path(imgdir, paste("corr_", regressand, ".pdf", sep = "")), width = 4 / 2.54, height = 4 / 2.54)
  # Create the heatmap
  # Create the heatmap
  p <- ggplot(cor_melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)),
      color = "black",
      size = 1.3,
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson's r",
      breaks = seq(-1, 1, 0.5)
    ) +
    guides(
      size = "none",
      colour = guide_colourbar(title.position = "right")
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev) +
    coord_fixed() +
    theme_minimal(base_size = 6) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 0, vjust = 0.5,
        color = "black", size = 3.5, face = "bold"
      ),
      axis.text.y = element_text(color = "black", size = 3.5, face = "bold"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.title = element_text(size = 3, face = "bold", angle = 90),
      legend.title.position = "left",
      legend.title.align = 0.5,
      legend.direction = "vertical",
      legend.text = element_text(size = 3),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.1, "cm"),
      plot.margin = margin(0.1, 0.1, 0.1, 0.1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  # Display the plot
  print(p)
  dev.off()
}



combined_estimate_plot <- function(model1, model2, filename, with_intercept = TRUE) {
  fixed1 <- extract_fixed(model1, with_intercept = with_intercept)
  rownames(fixed1) <- c("TDI", "RTT length (RBD)", "RTT length (NTD)")
  fixed2 <- extract_fixed(model2, with_intercept = with_intercept)
  rownames(fixed2) <- c("TDI", "RTT length (RBD)", "RTT length (NTD)")
  combined_fixed <- bind_rows(
    fixed1 %>%
      rownames_to_column(var = "Estimand") %>%
      mutate(dataset = "aa"),
    fixed2 %>%
      rownames_to_column(var = "Estimand") %>%
      mutate(dataset = "3di")
  )
  combined_fixed$term <- combined_fixed$Estimand
  pdf(filename,
    width = 4.63,
    height = 4.63
  )
  print(
    ggplot(combined_fixed, aes(y = term, x = Estimate, color = term, alpha = dataset)) +
      geom_point(size = 3, position = position_dodge(width = 0.25)) +
      geom_errorbar(aes(xmin = l95, xmax = u95),
        width = 0.2,
        position = position_dodge(width = 0.25)
      ) +
      labs(
        y = "",
        x = "Estimate"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_publish() +
      scale_color_ucscgb() +
      scale_alpha_manual(values = c("aa" = 1.0, "3di" = 0.5), name = "Dataset") +
      guides(color = "none") +
      theme(
        legend.position = "top",
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")
      )
  )
  dev.off()
}


combined_corr_plot <- function(cormat1, cormat2, outfile) {
  cols <- c("RTT length (RBD)", "RTT length (NTD)", "TDI")
  # Create first correlation matrix
  cor_matrix1 <- matrix(cormat1,
    nrow = 3,
    byrow = TRUE,
    dimnames = list(cols, cols)
  )

  # Create second correlation matrix
  cor_matrix2 <- matrix(cormat2,
    nrow = 3,
    byrow = TRUE,
    dimnames = list(cols, cols)
  )

  # Keep lower triangle from matrix 1, upper triangle from matrix 2
  cor_matrix1[upper.tri(cor_matrix1)] <- NA
  cor_matrix2[lower.tri(cor_matrix2, diag = TRUE)] <- NA

  # Convert to long format
  cor_melt1 <- melt(cor_matrix1)
  cor_melt2 <- melt(cor_matrix2)
  colnames(cor_melt1) <- c("Var1", "Var2", "value")
  colnames(cor_melt2) <- c("Var1", "Var2", "value")

  # Remove NA values and add identifier
  cor_melt1 <- cor_melt1[!is.na(cor_melt1$value), ]
  cor_melt1$matrix <- "lower"

  cor_melt2 <- cor_melt2[!is.na(cor_melt2$value), ]
  cor_melt2$matrix <- "upper"

  # Combine both datasets
  cor_melt <- rbind(cor_melt1, cor_melt2)

  pdf(outfile, width = 5 / 2.54, height = 4 / 2.54)

  # Create the split heatmap with two color scales
  p <- ggplot(cor_melt, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(
      data = subset(cor_melt, matrix == "lower"),
      aes(fill = value), color = "white", linewidth = 0.5
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "AA",
      breaks = seq(-1, 1, 0.5),
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        order = 1
      )
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = subset(cor_melt, matrix == "upper"),
      aes(fill = value), color = "white", linewidth = 0.5
    ) +
    scale_fill_gradient2(
      low = "#7B3294",
      mid = "white",
      high = "#008837",
      midpoint = 0,
      limits = c(-1, 1),
      name = "3Di",
      breaks = seq(-1, 1, 0.5),
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        order = 2
      )
    ) +
    geom_text(aes(label = sprintf("%.2f", value)),
      color = "black",
      size = 1.3,
      fontface = "bold"
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev) +
    coord_fixed() +
    theme_minimal(base_size = 6) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 0, vjust = 0.5,
        color = "black", size = 3.5, face = "bold"
      ),
      axis.text.y = element_text(color = "black", size = 3.5, face = "bold"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "right",
      legend.box = "horizontal",
      legend.title = element_text(size = 3, face = "bold"),
      legend.text = element_text(size = 3),
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.1, "cm"),
      plot.margin = margin(0.1, 0.1, 0.1, 0.1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  print(p)
  dev.off()
}
