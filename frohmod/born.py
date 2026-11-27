"""The BORN file (phonopy format) and the crystal symmetry needed to use it."""
import numpy as np
import spglib

from .supercell import find


def symmetry(cell):
    """Symmetry operations (fractional) and, for each atom, its independent partner.
    spglib with phonopy's tolerance, so the independent atoms are the ones phonopy expects."""
    s = cell["symbols"]
    data = spglib.get_symmetry_dataset((cell["lattice"], cell["positions"], [s.index(x) for x in s]), symprec=1e-5)
    if hasattr(data, "rotations"):  # spglib >= 2.5 returns an object, older versions a dict
        return data.rotations, data.translations, [int(j) for j in data.equivalent_atoms]
    return data["rotations"], data["translations"], [int(j) for j in data["equivalent_atoms"]]


def moved(cell, rot, tra, j):
    """Index of the atom that atom j lands on under x -> R x + t."""
    return find(rot @ cell["positions"][j] + tra, cell["symbols"][j], cell["lattice"], cell["positions"], cell["symbols"])


def cart(cell, rot):
    """A fractional rotation written in Cartesian coordinates."""
    l = cell["lattice"].T
    return l @ rot @ np.linalg.inv(l)


def read_born(path, cell):
    """BORN lists eps_inf and Z* of the independent atoms only; the others get the rotated Z*:
    Z_i = R Z_j R^T, with R the operation that carries atom j onto atom i."""
    with open(path) as f:
        lines = [x.split() for x in f if x.strip()]
    try:
        eps = np.array(lines[1][:9], float).reshape(3, 3)
        z = np.array([x[:9] for x in lines[2:]], float).reshape(-1, 3, 3)
    except (IndexError, ValueError):
        raise ValueError(f"{path} is not a BORN file (line 2: eps_inf, then one line of 9 numbers per Z*)")

    # check that the file fits this structure (and the atoms named on line 1, if it names them)
    rots, tras, eq = symmetry(cell)
    ind = sorted(set(eq))
    need = [j + 1 for j in ind]
    named = []
    if lines[0][0] == "#":
        named = [int(x) for x in lines[0] if x.isdigit()]
    if len(z) != len(ind) or (named and named != need):
        raise ValueError(f"{path} does not fit this structure: it needs one Z* line for each of atoms {need}")

    zj = dict(zip(ind, z))
    born = np.zeros((len(eq), 3, 3))
    for i, j in enumerate(eq):
        for rot, tra in zip(rots, tras):
            if moved(cell, rot, tra, j) == i:
                r = cart(cell, rot)
                born[i] = r @ zj[j] @ r.T
                break
    return eps, born


def write_born(path, cell, born, eps):
    """Like phonopy-vasp-born: make sum(Z*) = 0, average over symmetry, keep the independent atoms.
    born: Z* of every atom of `cell`; eps: eps_inf."""
    rots, tras, eq = symmetry(cell)
    rc = [cart(cell, r) for r in rots]
    born = born - born.mean(axis=0)
    eps = np.mean([r @ eps @ r.T for r in rc], axis=0)
    ind = sorted(set(eq))
    with open(path, "w") as f:
        f.write("# epsilon and Z* of atoms " + " ".join(str(j + 1) for j in ind) + "\n")
        f.write(" ".join(f"{x:12.8f}" for x in eps.ravel()) + "\n")
        for j in ind:
            # rotate the Z* of every atom equivalent to j back onto j, and average
            zs = [r.T @ born[moved(cell, rot, tra, j)] @ r for rot, tra, r in zip(rots, tras, rc)]
            f.write(" ".join(f"{x:12.8f}" for x in np.mean(zs, axis=0).ravel()) + "\n")
