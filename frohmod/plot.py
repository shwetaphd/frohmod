"""The two figures, drawn from modes.dat files only (one file = one phase, all its directions).

    plot 1, modes.png:    one panel per phase, a bar at every LO mode, height |g|^2 x n (meV), split by element
    plot 2, elements.png: one bar per phase, each element's share (%) over all LO modes of all directions

frohmod --plot draws them for its own run. To compare phases, run frohmod in each phase's folder, then:

    python -m frohmod.plot BiOBr/modes.dat Bi4O5Br2/modes.dat Bi3O4Br/modes.dat Bi5O7Br/modes.dat

The label of each phase is the name of the folder its modes.dat is in.
"""
import os
import sys

import numpy as np

# colours of the elements; other elements get matplotlib's default colours
colors = {"Bi": "darkturquoise", "Br": "green", "O": "red"}


def read_dat(path):
    """Label (folder name), elements, LO frequencies, |g|^2 x n and each element's part (meV) from a modes.dat."""
    with open(path) as f:
        names = f.readlines()[1].lstrip("# ").split()
    data = np.loadtxt(path, ndmin=2)
    col = {n: data[:, i] for i, n in enumerate(names)}
    els = [n[5:-4] for n in names if n.startswith("part_")]
    label = os.path.basename(os.path.dirname(os.path.abspath(path)))
    return {"label": label, "elements": els, "f_lo": col["LO_THz"], "g2n": col["g2n_meV"],
            "part": {el: col[f"part_{el}_meV"] for el in els}}


def share(d):
    """Each element's share (%) of |g|^2 x n over all LO modes in the file."""
    return {el: 100 * d["part"][el].sum() / d["g2n"].sum() for el in d["elements"]}


def color(el, k):
    return colors.get(el, f"C{k}")


def plot_modes(data, path="modes.png", temperature=300):
    """Plot 1: one panel per phase, every LO mode of every direction, |g|^2 x n split by element."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    els = data[0]["elements"]
    fig, axes = plt.subplots(len(data), 1, figsize=(6, 1 + 2 * len(data)), sharex=True, squeeze=False)
    for ax, d in zip(axes[:, 0], data):
        bottom = np.zeros(len(d["f_lo"]))
        for k, el in enumerate(els):
            ax.bar(d["f_lo"], d["part"][el], width=0.15, bottom=bottom, color=color(el, k))
            bottom += d["part"][el]
        ax.text(0.98, 0.9, d["label"], transform=ax.transAxes, ha="right", va="top", fontsize=9)
        ax.set(ylim=(0, 1.25 * bottom.max()))
    axes[-1, 0].set(xlim=(0, 1.05 * max(d["f_lo"].max() for d in data)), xlabel="LO phonon frequency (THz)")
    axes[0, 0].legend([plt.Rectangle((0, 0), 1, 1, color=color(el, k)) for k, el in enumerate(els)], els,
                      loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=len(els), frameon=False)
    fig.supylabel(f"$|g_\\nu|^2 \\times n_{{LO}}$ at {temperature:.0f} K (meV)")
    fig.savefig(path, dpi=300, bbox_inches="tight")


def plot_elements(data, path="elements.png", temperature=300):
    """Plot 2: one bar per phase, each element's share (%) of |g|^2 x n."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    els = data[0]["elements"]
    shares = [share(d) for d in data]
    fig, ax = plt.subplots(figsize=(6, 1 + 0.7 * len(data)))
    y, left = np.arange(len(data)), np.zeros(len(data))
    for k, el in enumerate(els):
        vals = np.array([s.get(el, 0) for s in shares])
        ax.barh(y, vals, left=left, color=color(el, k), label=el)
        for j in y:
            if vals[j] > 6:
                ax.text(left[j] + vals[j] / 2, j, f"{vals[j]:.0f}%", ha="center", va="center", color="white", fontsize=11)
        left += vals
    ax.set_yticks(y, [d["label"] for d in data])
    ax.invert_yaxis()
    ax.set(xlim=(0, 100), xlabel=f"Contribution in Fröhlich coupling at {temperature:.0f} K (%)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=len(els), frameon=False)
    fig.savefig(path, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    files = sys.argv[1:]
    if not files:
        sys.exit(__doc__)
    data = [read_dat(f) for f in files]
    plot_modes(data)
    plot_elements(data)
    print("saved modes.png (plot 1) and elements.png (plot 2)")
