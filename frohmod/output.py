"""Results as a printed table, a .dat file for plotting, and V_Sim animations (figures: plot.py)."""
import numpy as np


def table(r, eps, eps_ion=None, cutoff=3.0):
    """The table of one direction printed on screen. Modes below `cutoff` % of that direction's
    |g|^2 x n are not listed."""
    els, q = r["elements"], r["q"]
    out = [f"q = ({q[0]:.2f}, {q[1]:.2f}, {q[2]:.2f})   T = {r['temperature']:.0f} K   "
           f"eps_inf(q) = {r['eps_q']:.2f}   eps_0(q) = {r['eps0_q']:.2f}"]

    # checks: eps_0 against VASP, Z* sum, imaginary modes
    if eps_ion is not None:
        e0 = q @ (eps + eps_ion) @ q
        out.append(f"cross-check: eps_0(q) from VASP's epsilon_ion = {e0:.2f}")
        if abs(r["eps0_q"] - e0) > 0.1 * e0:
            out.append("warning: eps_0 differs from VASP's by > 10 %: BORN may not belong to this vasprun.xml")
    if r["born_sum"] > 0.1:
        out.append(f"warning: Z* add up to {r['born_sum']:.2f} e instead of 0: check the BORN file")
    if r["n_imaginary"]:
        out.append(f"warning: {r['n_imaginary']} imaginary mode(s), left out")

    # one row per LO mode, largest |g|^2 x n first
    head = f"  {'LO (THz)':>8} <- {'Gamma (THz)':>11}  {'mode':>4}  {'|g|^2 (meV)':>11}  {'n_LO':>5}  {'|g|^2 x n (meV)':>15}"
    head += "".join(f"  {el + ' (meV)':>8}" for el in els)
    out += ["", head]
    total = r["g2n"].sum()
    for v in np.argsort(-r["g2n"]):
        if r["g2n"][v] < cutoff / 100 * total:
            break
        k = r["gamma_mode"][v]
        line = f"  {r['f_lo'][v]:8.2f} <- {r['f_gamma'][k]:11.2f}  {k + 1:4d}  {r['g2'][v]:11.2f}  {r['n'][v]:5.2f}  "
        line += f"{r['g2n'][v]:15.2f}" + "".join(f"  {r['part'][el][v]:8.2f}" for el in els)
        out.append(line)
    out.append(f"  {'sum over all LO modes':>49}  {total:15.2f}" + "".join(f"  {r['part'][el].sum():8.2f}" for el in els))
    return "\n".join(out)


notes = """Gamma (THz): the same vibration at Gamma without the long-range field (no NAC)
mode: its phonopy mode number at Gamma (use it with --anime)
|g|^2: Frohlich coupling per phonon, (hbar w / 2) x Delta
n_LO: phonons in the mode at T
Bi, Br, ... (meV): each element's part of |g|^2 x n (|g|^2 x n x its share of the mode's vibration)"""


def save_dat(results, path):
    """Every mode of every direction in plain columns, ready to plot. part_X = |g|^2 x n x motion_X
    are the stacked bar pieces (meV); element X's share = sum of part_X / sum of g2n_meV."""
    els = results[0]["elements"]
    names = ["qx", "qy", "qz", "LO_THz", "Gamma_THz", "mode", "g2_meV", "n_LO", "g2n_meV"]
    names += [f"motion_{el}" for el in els] + [f"part_{el}_meV" for el in els]
    rows = []
    for r in results:
        cols = [np.tile(r["q"], (len(r["f_lo"]), 1)), r["f_lo"], r["f_gamma"][r["gamma_mode"]], r["gamma_mode"] + 1,
                r["g2"], r["n"], r["g2n"]]
        cols += [r["motion"][el] for el in els] + [r["part"][el] for el in els]
        rows.append(np.column_stack(cols))
    header = f"frohmod  T = {results[0]['temperature']:.0f} K  (one block of modes per direction q)\n" + " ".join(names)
    np.savetxt(path, np.vstack(rows), fmt=["%7.3f"] * 3 + ["%10.4f", "%10.4f", "%8d"] + ["%12.5f"] * (len(names) - 6),
               header=header)


def write_anime(path, cell, r, mode):
    """V_Sim file in phonopy's anime.ascii format for mode `mode` (phonopy's numbering at Gamma, without NAC).
    It holds what that mode becomes along q, i.e. the LO mode with the largest overlap, as the
    displacement u = e / sqrt(M) (phonopy's own anime at Gamma only gives the mode without NAC)."""
    # V_Sim wants the cell turned so that a is along x and b in the xy plane;
    # the rows of the Cholesky factor of the metric L L^T are exactly those vectors
    lat = cell["lattice"]
    a = np.linalg.cholesky(lat @ lat.T)

    t = mode - 1
    v = int(r["overlap"][t].argmax())  # the LO mode this mode turns into along q
    u = r["e_lo"][:, v].reshape(-1, 3) / np.sqrt(cell["masses"])[:, None] @ np.linalg.inv(lat) @ a
    q = r["q"]
    lines = [f"# frohmod for v_sim: mode {mode} ({r['f_gamma'][t]:.4f} THz) -> {r['f_lo'][v]:.4f} THz "
             f"along q = ({q[0]:.3f}, {q[1]:.3f}, {q[2]:.3f})",
             f"{a[0, 0]:15.9f}{a[1, 0]:15.9f}{a[1, 1]:15.9f}",
             f"{a[2, 0]:15.9f}{a[2, 1]:15.9f}{a[2, 2]:15.9f}"]
    for (x, y, z), s in zip(cell["positions"] @ a, cell["symbols"]):
        lines.append(f"{x:15.9f}{y:15.9f}{z:15.9f} {s:>2s}")
    lines.append(f"#metaData: qpt=[0.000000;0.000000;0.000000;{r['f_lo'][v]:f} \\")
    for d in u:
        lines.append("#; " + "; ".join(f"{x:f}" for x in list(d) + [0.0, 0.0, 0.0]) + " \\")
    lines.append("# ]")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
