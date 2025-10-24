import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from src.iqtreestats import get_tree_stats, get_tree_weights


def clear_axes(
    ax=None, top=True, right=True, left=False, bottom=False, minorticks_off=True
):
    """
    A more forcing version of sns.despine.

    From: https://gist.github.com/Neclow/91867667fdc463c0f36aaa7a26e84e0e

    ax : matplotlib axes, optional
        Specific axes object to despine. Ignored if fig is provided.
    top, right, left, bottom : boolean, optional
        If True, remove that spine.
    minorticks_off: boolean, optional
        If True, remove all minor ticks
    """
    if ax is None:
        axes = plt.gcf().axes
    else:
        axes = [ax]

    for ax_i in axes:
        sns.despine(ax=ax_i, top=top, right=right, left=left, bottom=bottom)
        if minorticks_off:
            ax_i.minorticks_off()
        ax_i.tick_params(axis="x", which="both", top=not top)
        ax_i.tick_params(axis="y", which="both", right=not right)
        ax_i.tick_params(axis="y", which="both", left=not left)
        ax_i.tick_params(axis="x", which="both", bottom=not bottom)


def set_size(width, layout="h", fraction=1):
    """Set figure dimensions to avoid scaling in LaTeX.

    Heavily inspired by: https://jwalton.info/Embed-Publication-Matplotlib-Latex/
    Parameters
    ----------
    width: float
        Document textwidth or columnwidth in pts
        Report: 390 pt
    layout: string
        h: horizontal layout
        v: vertical layout
        s: square layout
    fraction: float, optional
        Fraction of the width which you wish the figure to occupy
    Returns
    -------
    fig_dim: tuple
        Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio for aesthetic figures
    # https://disq.us/p/2940ij3
    golden_ratio = (5**0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    if layout == "h":
        fig_height_in = fig_width_in * golden_ratio
    elif layout == "v":
        fig_height_in = fig_width_in / golden_ratio
    elif layout == "s":
        fig_height_in = fig_width_in
    else:
        raise ValueError(f"Unknown layout: {layout}")

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim


def plot_tree_stats(
    tr_path, max_values, *, outgroup=None, title=None, padding=1.0, order=None
):
    df = get_tree_stats(tr_path, outgroup=outgroup)
    df["weight"] = get_tree_weights(tr_path.replace(".treefile.annotated", ".iqtree"))
    df = df.loc[:, list(max_values.keys())].reset_index()
    print(df.head())
    if order is not None:
        df = df.iloc[order, :]
    print(df.head())

    id_column = "index"
    categories = df._get_numeric_data().columns.tolist()
    data = df[categories].to_dict(orient="list")
    ids = df[id_column].tolist()

    max_values = {key: padding * value for key, value in max_values.items()}
    normalized_data = {
        key: np.array(values) / max_values[key] for key, values in data.items()
    }
    n_vars = len(data.keys())
    tiks = list(data.keys())
    tiks += tiks[:1]
    angles = np.linspace(np.pi / 4, 2 * np.pi, n_vars, endpoint=False).tolist() + [
        np.pi / 4
    ]
    angles = list(np.array([45, 135, 225, 315]) * np.pi / 180) + [np.pi / 4]
    fig, ax = plt.subplots(figsize=set_size(260, "s"), subplot_kw=dict(polar=True))
    texts = []
    palette = sns.color_palette("colorblind").as_hex()

    # Plot the data
    for i, model_name in enumerate(ids):
        values = [normalized_data[key][i] for key in data.keys()]
        actual_values = [data[key][i] for key in data.keys()]
        values += values[:1]  # Close the plot for a better look
        ax.plot(angles, values, label=model_name, color=palette[i])

    ax.legend(title="Class", loc="upper left", bbox_to_anchor=(1.1, 1.0), frameon=True)

    # Add the central 0 tick manually
    text = ax.text(
        0,
        0,
        "0",
        ha="center",
        va="center",
        fontsize=10,
        color="gray",
    )
    texts.append(text)
    # adjust_text(texts)

    for i, cat in enumerate(tiks[:-1]):
        angle = angles[i]

        yticks = ax.get_yticks()
        tick_positions = [tick for tick in yticks if 0 < tick <= 1]

        # Add labels for key tick positions
        for tick_pos in tick_positions[::2]:  # Every other tick to avoid crowding
            actual_value = tick_pos * max_values[cat]

            # Calculate slight offset angle to avoid overlap with axis line
            offset_angle = angle + 0.01  # Smaller angular offset

            text = ax.text(
                offset_angle,
                tick_pos,
                f"{actual_value:.2f}",
                ha="center",
                va="top" if i < len(tiks) - 1 else "bottom",
                fontsize=9,  # Smaller font
                color="gray",
            )
            texts.append(text)

    # Clean up the axes
    ax.set_ylim(0, 1)  # Force the radial axis to go from 0 to 1
    ax.set_yticklabels([])  # Hide default radial tick labels
    ax.set_xticks(angles[:-1])  # Remove duplicate angle
    ax.set_xticklabels(tiks[:-1])  # Remove duplicate label
    labels = [tick.get_text() for tick in ax.get_xticklabels()]
    pretty_labels = [
        lbl if lbl.isupper() else lbl.capitalize().replace("_", "\n") for lbl in labels
    ]
    ax.set_xticklabels(pretty_labels)
    ax.tick_params(axis="x", which="major", pad=15)
    ax.yaxis.grid(alpha=0.2, zorder=0)
    ax.xaxis.grid(alpha=0.2, zorder=0)  # Adds subtle radial lines
    ax.spines["polar"].set_visible(False)

    if title is not None:
        plt.suptitle(title)
    return fig, ax


def plot_mastscores(mast_scores):
    """
    Plot scores of MAST runs for different classes and criteria.

    Parameters
    ----------
    mast_scores : pd.DataFrame
        Columns = ["n_classes", "AIC", "AICc", "BIC"]
        Rows = Scores for different number of classes.
    """
    mast_melt = (
        mast_scores.reset_index()
        .melt(id_vars=["n_classes"], var_name="Criterion", value_name="Score")
        .sort_values(by="n_classes")
    )

    palette = sns.color_palette("Blues_r", n_colors=mast_scores.shape[0])
    if len(mast_melt.Criterion.unique()) == 1:
        fig, ax = plt.subplots(figsize=set_size(260, "h"))
        sns.scatterplot(
            x="n_classes",
            y="Score",
            data=mast_melt,
            s=50,
            hue="Score",
            palette=palette,
            edgecolor="k",
        )
        ax.set_ylabel(mast_melt.Criterion.unique()[0])
        ax.get_legend().remove()
    else:
        fig, ax = plt.subplots(figsize=set_size(280, "h"))
        sns.stripplot(
            x="Criterion",
            hue="n_classes",
            order=sorted(mast_melt["Criterion"].unique()),
            y="Score",
            dodge=True,
            size=5,
            data=mast_melt,
            palette=palette,
        )
        plt.legend(
            bbox_to_anchor=(1.01, 1.0),
            loc="upper left",
            title=r"\# Classes",
            fontsize=10,
            title_fontsize=11,
            alignment="left",
        )
    ax.grid(alpha=0.2)

    ax.set_xlabel(r"\# Classes")
    plt.tight_layout()
    clear_axes()
    return fig, ax


def plot_sites(slh, domains, palette, figname, label, offset=0, order=None):
    # slh.plot()
    # plt.show()

    w, h = set_size(290, "h")
    fig4, ax4 = plt.subplots(1, 1, figsize=(2 * w, h))

    # Get the data without the dropped column
    if r"$L_{MAST}$" in slh.columns:
        # data = slh.sub(slh.loc[:, r"$L_{MAST}$"], axis=0).drop(r"$L_{MAST}$", axis=1, errors="ignore")
        # Subtract mean from each row
        data = slh.drop(r"$L_{MAST}$", axis=1, errors="ignore")
        # data = data.sub(slh.mean(axis=1), axis=0)

    ymax = offset

    # Sort columns if order is not None
    if order is not None:
        data = data.loc[:, order]

    # Don't use the built-in plot method, instead plot manually with varying alpha
    x_data = data.index

    # For each column (line), plot segments with varying alpha
    for col_idx, (col_name, col_data) in enumerate(data.items()):
        x_vals = []
        y_vals = []
        alphas = []

        for i in range(len(x_data)):
            # Get y values for all columns at this x position
            y_values_at_x = [data.iloc[i, j] for j in range(len(data.columns))]
            max_idx = np.argmax(y_values_at_x)

            # Set alpha based on whether this line has the maximum value at this x
            if col_idx == max_idx:
                x_vals.append(x_data[i])
                other_vals = y_values_at_x[:col_idx] + y_values_at_x[col_idx + 1 :]
                y_vals.append(col_data.iloc[i] - np.mean(other_vals))
                alphas.append(1.0)  # Opaque for maximum
            else:
                x_vals.append(x_data[i])
                y_vals.append(0)
                alphas.append(0.01)  # Transparent for non-maximum

        # Plot line segments with varying alpha
        for i in range(len(x_vals) - 1):
            ax4.plot(
                [x_vals[i], x_vals[i + 1]],
                [y_vals[i], y_vals[i + 1]],
                color=palette[col_idx],
                alpha=min(
                    alphas[i], alphas[i + 1]
                ),  # Use minimum alpha of the two points
                label=col_name if i == 0 else "",
            )  # Only label once

    # Add legend
    handles, labels = ax4.get_legend_handles_labels()
    for handle in handles:
        handle.set_alpha(1.0)  # Make legend lines very transparent
    ax4.legend(handles, labels, bbox_to_anchor=(0.97, 1), loc="upper left")

    clear_axes()
    ax4.set_ylabel(label)

    for i, dom in enumerate(["S1 (NTD)", "S2", "S1 (RBD)"]):
        start, end = domains[dom]
        ax4.axvspan(
            start,
            end,
            color="lightgrey",
            alpha=0.3,
            ymax=1.0,
            clip_on=False,
            zorder=0,
        )
        ax4.text(
            start + (end - start) / 2,
            ymax,
            dom,
            weight="semibold",
            fontsize=7,
            ha="left",
            va="bottom",
            rotation=45,
        )
    plt.savefig(f"img/{figname}.pdf", bbox_inches="tight", pad_inches=0.1)
    plt.show()
