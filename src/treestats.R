library(ape)
library(phangorn)
library(phytools)

pathnode <- function(phylo, tipsonly = TRUE) {
  dist_nodes <- ape::dist.nodes(phylo)
  root <- length(phylo$tip.label) + 1 # phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]

  if (tipsonly == TRUE) {
    roottotippath <- dist_nodes[
      as.numeric(rownames(dist_nodes)) == root,
      seq_along(phylo$tip.label)
    ]
    nodesinpath <- sapply(seq_along(phylo$tip.label), function(x) {
      length(Ancestors(
        phylo,
        x
      ))
    })
  } else {
    roottotippath <- dist_nodes[as.numeric(rownames(dist_nodes)) == root, ]
    nodesinpath <- sapply(
      1:(length(phylo$tip.label) + phylo$Nnode),
      function(x) length(Ancestors(phylo, x))
    )
  }

  list(
    roottotippath = roottotippath,
    nodesinpath = nodesinpath
  )
}

rtt_cov <- function(tre) {
  if (!is.rooted(tre)) tre <- midpoint(tre)
  rtts <- pathnode(tre)[[1]]
  rttcov <- sd(rtts, na.rm = TRUE) / mean(rtts, na.rm = TRUE)
  rttcov
}

treestats.single <- function(tr, outgroup = NULL) {
  if (!is.null(outgroup)) {
    tr <- root(tr, outgroup)
  }
  if (inherits(tr, "multiPhylo") && length(tr) > 1) {
    stop("Multiple trees found, please provide a single tree.")
  }
  if (is.null(tr$node.label)) {
    meanshalrt <- NA
    meanufboot <- NA
  } else {
    split_labels <- strsplit(tr$node.label, "/")
    nsplit <- length(split_labels[[2]])
    if (nsplit == 0) {
      meanshalrt <- NA
      meanufboot <- NA
    } else if (nsplit == 1) {
      meanshalrt <- NA
      meanufboot <- mean(
        sapply(split_labels, function(x) as.numeric(x[1])),
        na.rm = TRUE
      )
    } else if (nsplit == 2) {
      meanshalrt <- mean(
        sapply(split_labels, function(x) as.numeric(x[1])),
        na.rm = TRUE
      )
      meanufboot <- mean(
        sapply(split_labels, function(x) as.numeric(x[2])),
        na.rm = TRUE
      )
    } else {
      stop("Unexpected number of splits in node labels.")
    }
  }
  leaf_idxs <- which(tr$edge[, 2] <= Ntip(tr))
  tip_length <- sum(tr$edge.length[leaf_idxs], na.rm = TRUE)
  total_length <- sum(tr$edge.length, na.rm = TRUE)
  internal_length <- total_length - tip_length

  data.frame(
    meanshalrt = meanshalrt,
    meanufboot = meanufboot,
    MCL = rtt_cov(tr),
    tree_length = total_length,
    internal_length = internal_length,
    tip_length = tip_length,
    stemminess = internal_length / total_length,
    tips = Ntip(tr)
  )
}

treestats <- function(file, outgroup = NULL) {
  tr <- read.tree(file)
  if (inherits(tr, "multiPhylo") && length(tr) > 1) {
    stats_list <- lapply(tr, function(single_tr) treestats.single(single_tr, outgroup))
    out <- do.call(rbind, stats_list)
  } else {
    out <- treestats.single(tr, outgroup)
  }
  out
}
