import subprocess

import numpy as np
import pandas as pd
import rpy2.robjects as ro

from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from scipy.ndimage import uniform_filter1d


def get_tree_weights(iqtree_file):
    command = f"grep 'Tree weights' {iqtree_file} | cut -d ':' -f2"
    result = subprocess.run(
        command, shell=True, capture_output=True, text=True, check=True
    )
    return list(map(float, result.stdout.strip().split(",")))


def get_tree_stats(tr_path, outgroup=None):
    """Extract shape statistics from a phylogenetic tree"""
    with localconverter(ro.default_converter + pandas2ri.converter):
        importr("phangorn")
        ro.globalenv["file"] = tr_path
        ro.globalenv["outgroup"] = outgroup or ro.NULL
        tr_stats = ro.r(
            """
            source("src/treestats.R")
            treestats(file, outgroup = outgroup)
            """
        ).replace(np.iinfo(np.int32).min, float("nan"))

    return tr_stats


def prepare_sitelh(sitelh_file, window):
    slh = pd.read_csv(sitelh_file, index_col=0, header=0, sep="\t", comment="#")

    slh.columns = [r"$L_{MAST}$", *[rf"$L_{i}$" for i in range(1, slh.shape[1])]]

    arr = slh.to_numpy()

    arr = uniform_filter1d(arr, size=window, axis=0, mode="nearest")

    return pd.DataFrame(arr, index=slh.index, columns=slh.columns)
