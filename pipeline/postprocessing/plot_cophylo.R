library(ape)
library(dplyr)
library(ggsci)
library(phytools)
library(Cairo)

source("pipeline/plot_config.R")

create_color_palette <- function(virusinfo) {
  subgenera <- unique(virusinfo$Subgenus)
  palette <- pal_npg("nrc")(length(subgenera))
  names(palette) <- subgenera
  palette
}

create_color_mapping <- function(tr, virusinfo, palette) {
  colors <- c()
  for (i in seq_along(tr$tip.label)) {
    label <- tr$tip.label[[i]]
    parts <- strsplit(label, "/")[[1]]
    virus_name <- paste(head(parts, -2), collapse = "/")
    virus_name <- gsub("'", "", virus_name)
    virus_name <- gsub(" ", "", virus_name)
    virus_name <- gsub("Cov", "CoV", virus_name)
    subgenus <- virusinfo %>%
      filter(Virus == virus_name) %>%
      select(Subgenus) %>%
      slice(1) %>%
      pull()
    colors <- c(colors, palette[[subgenus]])
  }
  colors
}

create_shape_palette <- function(virusinfo) {
  receptors <- sort(unique(virusinfo$Receptor))
  symbols <- c(
    "\u25CF", "\u25A0", "\u25C6", "\u25B2",
    "\u25A1", "\u25B3", "\u25C7", "\u25CB"
  )
  palette <- head(symbols, length(receptors))
  names(palette) <- receptors
  palette
}

create_shape_mapping <- function(tr, virusinfo, palette, show_tip_label = TRUE, right = TRUE) {
  shapes <- c()
  for (i in seq_along(tr$tip.label)) {
    label <- tr$tip.label[[i]]
    parts <- strsplit(label, "/")[[1]]
    virus_name <- paste(head(parts, -2), collapse = "/")
    # Remove ' and spaces if present
    virus_name <- gsub("'", "", virus_name)
    virus_name <- gsub(" ", "", virus_name)
    virus_name <- gsub("Cov", "CoV", virus_name)
    receptor <- virusinfo %>%
      filter(Virus == virus_name) %>%
      select(Receptor) %>%
      slice(1) %>%
      pull()
    new_label <- if (show_tip_label && right) {
      paste0(label, "  ", palette[[receptor]])
    } else if (show_tip_label && !right) {
      paste0(palette[[receptor]], "  ", label)
    } else {
      palette[[receptor]]
    }
    shapes <- c(shapes, new_label)
  }
  shapes
}

cophylo_plot <- function(tree1, tree2, outgroup, filename, show_tip_label = TRUE) {
  tree1 <- ladderize(root(tree1, outgroup = outgroup, resolve.root = TRUE))
  tree2 <- ladderize(root(tree2, outgroup = outgroup, resolve.root = TRUE))

  shape_palette <- create_shape_palette(virusinfo)
  if (show_tip_label) {
    tree1$tip.label <- create_shape_mapping(
      tree1, virusinfo,
      shape_palette, show_tip_label
    )
    tree2$tip.label <- create_shape_mapping(
      tree2, virusinfo,
      shape_palette, show_tip_label
    )
  }

  cophy <- cophylo(
    tree1, tree2,
    assoc = cbind(tree1$tip.label, tree1$tip.label),
    rotate = TRUE,
    use.edge.length = TRUE
  )

  color_palette <- create_color_palette(virusinfo)
  colors <- create_color_mapping(tree1, virusinfo, color_palette)

  width <- if (show_tip_label) 12 else 6
  CairoPDF(filename, family = "URWGothic", width = width, height = 8)
  par(mar = c(7, 4, 4, 2) + 0.1, family = "URWGothic") # increase bottom margin
  ftype <- if (show_tip_label) "reg" else "off"
  plot(cophy,
    link.type = "curved",
    link.lwd = 3,
    link.lty = "solid",
    link.col = colors,
    node.lwd = 0,
    tip.col = NULL,
    cex = 0,
    fsize = 0.66,
    scale.bar = rep(1, 2),
    pts = FALSE,
    ftype = ftype
  )
  if (!show_tip_label) {
    cophy$trees[[1]]$tip.label <- create_shape_mapping(
      cophy$trees[[1]],
      virusinfo,
      shape_palette,
      show_tip_label
    )
    cophy$trees[[2]]$tip.label <- create_shape_mapping(
      cophy$trees[[2]],
      virusinfo,
      shape_palette,
      show_tip_label
    )
    tiplabels.cophylo(
      text = cophy$trees[[1]]$tip.label,
      bg = "none",
      frame = "none",
      which = "left",
      cex = 1
    )
    tiplabels.cophylo(
      text = cophy$trees[[2]]$tip.label,
      bg = "none",
      frame = "none",
      which = "right",
      cex = 1
    )
  }

  # Add legend
  max_rows <- max(length(color_palette), length(shape_palette)) # 5

  palette_pad <- c(color_palette, rep(NA, max_rows - length(color_palette)))
  palette_names_pad <- c(
    names(color_palette),
    rep("", max_rows - length(color_palette))
  )

  shape_names_pad <- names(shape_palette)
  shape_pad <- shape_palette

  # Combine legends for a 2-column grid
  legend(
    "bottom",
    legend = c(palette_names_pad, shape_names_pad),
    col = c(palette_pad, rep("black", length(shape_pad))),
    lwd = c(ifelse(is.na(palette_pad), NA, 2), rep(NA, length(shape_pad))),
    pch = c(rep(NA, length(palette_pad)), shape_pad),
    pt.bg = c(rep(NA, length(palette_pad)), rep("white", length(shape_pad))),
    ncol = 2,
    title = "Legend",
    cex = 0.6,
    bty = "o",
    bg = "white"
  )
  dev.off()
}

cophylo_plot(
  ntd_tree_aa, ntd_tree_3di,
  outgroup = outgroup,
  filename = file.path(imgdir, "cophylo_ntd.pdf")
)

cophylo_plot(
  rbd_tree_aa, rbd_tree_3di,
  outgroup = outgroup,
  filename = file.path(imgdir, "cophylo_rbd.pdf")
)

cophylo_plot(
  s2_tree_aa, s2_tree_3di,
  outgroup = outgroup,
  filename = file.path(imgdir, "cophylo_s2.pdf")
)

# RBD vs NTD (AA)
cophylo_plot(
  rbd_tree_aa, ntd_tree_aa,
  outgroup = outgroup,
  filename = file.path(imgdir, "cophylo_aa_rbd_vs_ntd.pdf")
)
cophylo_plot(
  rbd_tree_aa, ntd_tree_aa,
  outgroup = outgroup,
  filename = file.path(imgdir, "cophylo_aa_rbd_vs_ntd_notip.pdf"),
  show_tip_label = FALSE
)

# RBD vs NTD (3Di)
cophylo_plot(
  rbd_tree_3di, ntd_tree_3di,
  outgroup = outgroup,
  filename = file.path(imgdir, "cophylo_3di_rbd_vs_ntd.pdf")
)
cophylo_plot(
  rbd_tree_3di, ntd_tree_3di,
  outgroup = outgroup,
  filename = file.path(imgdir, "cophylo_3di_rbd_vs_ntd_notip.pdf"),
  show_tip_label = FALSE
)
