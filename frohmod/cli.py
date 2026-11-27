"""Command line:

    frohmod --phonon vasprun.xml [--dielectric-born vasprun.xml] [--dim 2 2 1] [--q 1 0 0 --q 0 0 1]
              [-T 300] [--born BORN] [--anime 6 13 [--vesta structure.vesta]] [--plot [modes] [elements]]

--phonon          VASP run with the phonons (IBRION = 5-8): the primitive cell or a supercell
--dielectric-born VASP run with Born charges and eps_inf (LCALCEPS or LEPSILON); default: the --phonon run
--dim             supercell size of the --phonon run, as phonopy's DIM; default: found from the structure
--q               direction of the phonon wave, Cartesian; give it several times for several directions
                  (default: x, y and z)
-T                temperature in K (default 300)
--born            use this BORN file instead of writing one from --dielectric-born
--anime           write anime-<mode>_q<q>.ascii (V_Sim) for these mode numbers of the table (phonopy's numbering);
                  each file holds what that mode becomes along q (its LO displacement)
--vesta          with --anime: also draw each mode as arrows in VESTA, in anime-<mode>_q<q>.vesta;
                  give the POSCAR of the same calculation opened in VESTA and saved as a .vesta file
--plot            save the figures (needs matplotlib): "modes" = modes.png, |g|^2 x n of every LO mode;
                  "elements" = elements.png, each element's share; --plot alone = both.
                  To compare phases in one figure: python -m frohmod.plot */modes.dat

Writes to the current folder: BORN (unless --born is given), modes.dat (every mode of every direction,
for plotting), anime-<mode>_q<q>.ascii (and .vesta) and the figures if asked.
"""
import argparse
import os
import sys

import numpy as np

from .analysis import analyse, element_share
from .born import read_born, symmetry, write_born
from .output import notes, save_dat, table, write_anime
from .read import read_vasprun
from .supercell import born_on_unit_cell, unit_cell
from .vesta import ascii_to_vesta


def main(argv=None):
    p = argparse.ArgumentParser(prog="frohmod",
                                description="Mode-resolved analysis of Frohlich scattering from VASP phonon calculations.",
                                epilog="Writes to the current folder: BORN, modes.dat, anime-<mode>_q<q>.ascii")
    p.add_argument("--phonon", required=True, metavar="VASPRUN",
                   help="vasprun.xml with the phonons (IBRION = 5-8): primitive cell or supercell")
    p.add_argument("--dielectric-born", metavar="VASPRUN",
                   help="vasprun.xml with Born charges and eps_inf (LCALCEPS or LEPSILON); default: the --phonon file")
    p.add_argument("--dim", nargs="+", type=int, metavar="N",
                   help="supercell size of --phonon, as phonopy's DIM (3 or 9 integers); default: found from the structure")
    p.add_argument("--q", nargs=3, type=float, action="append", metavar=("X", "Y", "Z"),
                   help="direction of the phonon wave, Cartesian; repeat for more directions (default: x, y and z)")
    p.add_argument("-T", "--temperature", type=float, default=300.0, metavar="K", help="temperature in K (default 300)")
    p.add_argument("--born", metavar="FILE", help="use this BORN file (phonopy format) instead of writing one")
    p.add_argument("--anime", nargs="+", type=int, metavar="MODE",
                   help="write anime-<mode>.ascii (V_Sim) for these mode numbers of the table")
    p.add_argument("--vesta", metavar="FILE",
                   help="with --anime: also write anime-<mode>_q<q>.vesta (arrows), from this structure saved by VESTA")
    p.add_argument("--plot", nargs="*", choices=["modes", "elements"], metavar="FIGURE",
                   help="save modes.png (|g|^2 x n of every LO mode) and/or elements.png (each element's share); "
                        "--plot alone: both (needs matplotlib)")
    a = p.parse_args(argv)
    if a.q is None:
        a.q = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    if a.vesta and not a.anime:
        p.error("--vesta needs --anime (which modes to draw)")
    if a.dim is not None and len(a.dim) not in (3, 9):
        p.error("--dim takes 3 or 9 integers")

    try:
        # phonons: read the run and fold it onto the unit cell
        ph = read_vasprun(a.phonon)
        if ph["fc"] is None:
            raise ValueError(f"{a.phonon} has no phonons: give the vasprun.xml of an IBRION = 5-8 run")
        cell = unit_cell(ph, a.dim)
        if len(cell["symbols"]) != len(ph["symbols"]):
            print(f"supercell of {len(ph['symbols'])} atoms folded onto the {len(cell['symbols'])}-atom unit cell")
            if a.dim is None:
                print("--dim found from the structure: " + " ".join(map(str, cell["dim"].ravel())))
        n = sum(np.allclose(rot, np.eye(3)) for rot in symmetry(cell)[0])  # pure translations
        if n > 1:
            print(f"warning: this cell still holds {n} copies of a smaller cell: is --dim too small?")

        # Born charges and eps_inf: write BORN from the dielectric run, or use the one given
        src = ph
        if a.dielectric_born:
            src = read_vasprun(a.dielectric_born)
        if a.born:
            if not os.path.exists(a.born):
                raise ValueError(f"{a.born} not found")
            born_file = a.born
        else:
            if src["born"] is None or src["eps"] is None:
                raise ValueError(f"{src['path']} has no Born charges / dielectric tensor: give --dielectric-born "
                                 "(a run with LCALCEPS or LEPSILON), or --born")
            born_file = "BORN"
            write_born(born_file, cell, born_on_unit_cell(cell, src), src["eps"])
            print(f"wrote BORN from {src['path']}")

        # the analysis itself
        eps, born = read_born(born_file, cell)
        results = [analyse(cell, born, eps, q, a.temperature) for q in a.q]
        for m in a.anime or []:
            if not 1 <= m <= len(results[0]["f_lo"]):
                raise ValueError(f"--anime {m}: the modes are numbered 1 to {len(results[0]['f_lo'])}")
        if a.vesta and not os.path.exists(a.vesta):
            raise ValueError(f"{a.vesta} not found")
    except ValueError as e:
        sys.exit(f"error: {e}")

    # results: one table per direction, then each element's share over all of them
    for r in results:
        print("\n" + table(r, eps, src["eps_ion"]))
    share = element_share(results)
    qs = ", ".join("(" + " ".join(f"{x:g}" for x in q) + ")" for q in a.q)
    print("\n" + notes + "\n")
    print(f"Frohlich coupling at {a.temperature:.0f} K from all LO modes along {qs}:  "
          + "   ".join(f"{el} {share[el]:.1f} %" for el in share))

    # files: .dat, animations, figure
    save_dat(results, "modes.dat")
    print("\ndata for plotting: modes.dat")
    for r in results:
        tag = "_".join(f"{x:g}" for x in r["q"])
        for m in a.anime or []:
            name = f"anime-{m}_q{tag}"
            write_anime(name + ".ascii", cell, r, m)
            v = int(r["overlap"][m - 1].argmax())
            print(f"V_Sim animation: {name}.ascii "
                  f"(mode {m}, {r['f_gamma'][m - 1]:.2f} THz -> {r['f_lo'][v]:.2f} THz along q)")
            if a.vesta:
                try:
                    ascii_to_vesta(name + ".ascii", a.vesta, name + ".vesta")
                    print(f"VESTA arrows: {name}.vesta")
                except ValueError as e:
                    sys.exit(f"error: {e}")
    if a.plot is not None:
        from .plot import plot_elements, plot_modes, read_dat  # here, so python -m frohmod.plot runs cleanly
        figures = a.plot or ["modes", "elements"]
        data = [read_dat("modes.dat")]
        if "modes" in figures:
            plot_modes(data, "modes.png", a.temperature)
            print("figure: modes.png")
        if "elements" in figures:
            plot_elements(data, "elements.png", a.temperature)
            print("figure: elements.png")
