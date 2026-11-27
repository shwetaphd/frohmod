"""From the phonon run (primitive cell or a supercell) to the unit cell at Gamma.

Like phonopy's DIM: unit-cell lattice = supercell lattice / dim. If no dim is given, it is found from
the structure. Atoms are matched by position, so the order of atoms in the supercell does not matter.
"""
import itertools

import numpy as np

tol = 0.1  # A, for matching atoms between cells


def find(p, s, lattice, positions, symbols):
    """Index of the atom of type s at fractional position p (any lattice translation away), or None."""
    d = positions - p
    d -= np.round(d)
    r = np.linalg.norm(d @ lattice, axis=1)
    r[[x != s for x in symbols]] = np.inf
    k = int(np.argmin(r))
    if r[k] < tol:
        return k
    return None


def detect_dim(sc):
    """Supercell matrix m (3x3 integers) with supercell lattice = m @ unit-cell lattice,
    found from the structure itself (the same m that phonopy's DIM / --dim takes)."""
    lattice = np.asarray(sc["lattice"], float)
    pos = np.asarray(sc["positions"], float) % 1
    sym = sc["symbols"]

    # pure translations: shifts that map every atom onto an atom of the same type
    shifts = []
    for j, s in enumerate(sym):
        if s != sym[0]:
            continue
        t = pos[j] - pos[0]
        t -= np.round(t)
        if all(find(p + t, s2, lattice, pos, sym) is not None for p, s2 in zip(pos, sym)):
            shifts.append(t)
    n = len(shifts)  # number of unit cells in this cell
    if n == 1:
        return np.eye(3, dtype=int)

    # the translations form a lattice: take its three shortest independent vectors
    grid = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    cand = np.array([t + g for t in shifts for g in grid])
    cart = cand @ lattice
    length = np.linalg.norm(cart, axis=1)
    basis = []
    for i in np.argsort(length):
        if length[i] < 1e-6:
            continue
        trial = basis + [i]
        unit = cart[trial] / length[trial][:, None]
        if np.linalg.matrix_rank(unit, tol=1e-2) == len(trial):
            basis = trial
        if len(basis) == 3:
            break
    p = cand[basis]
    if len(basis) < 3 or abs(abs(np.linalg.det(p)) * n - 1) > 1e-2:
        raise ValueError("could not find the primitive cell of the --phonon cell: give --dim")
    if np.linalg.det(p) < 0:
        p[2] *= -1  # keep a right-handed cell
    m = np.linalg.inv(p)
    mi = np.rint(m).astype(int)
    if not np.allclose(m, mi, atol=1e-2):
        raise ValueError("the --phonon cell is not an integer supercell of its primitive cell: give --dim")
    return mi


def unit_cell(sc, dim=None):
    """Unit cell, with the force constants folded to Gamma: Phi(k, k') = (1/N) sum over all copies of k and k'.
    dim: 3 or 9 integers as phonopy's DIM, or None to find it from the structure."""
    if dim is None:
        m = detect_dim(sc)
    elif len(dim) == 3:
        m = np.diag(dim)
    else:
        m = np.reshape(dim, (3, 3))
    ncell = int(round(abs(np.linalg.det(m))))
    lattice = np.linalg.inv(m) @ sc["lattice"]
    frac = (sc["positions"] @ sc["lattice"] @ np.linalg.inv(lattice)) % 1

    # keep the first atom found at each unit-cell position; owner[i] = unit-cell atom that atom i copies
    pos, sym, mass, owner = [], [], [], []
    for p, s, w in zip(frac, sc["symbols"], sc["masses"]):
        k = None
        if pos:
            k = find(p, s, lattice, np.array(pos), sym)
        if k is None:
            k = len(pos)
            pos.append(p)
            sym.append(s)
            mass.append(w)
        owner.append(k)
    if len(pos) * ncell != len(owner):
        raise ValueError(f"the --phonon cell is not a {' '.join(map(str, dim or m.ravel()))} supercell: check --dim")

    # Phi_unit = C^T Phi_sc C / N, with C[i, k] = 1 when supercell atom i is a copy of unit-cell atom k
    c = np.zeros((3 * len(owner), 3 * len(pos)))
    for i, k in enumerate(owner):
        c[3*i:3*i + 3, 3*k:3*k + 3] = np.eye(3)
    return {"lattice": lattice, "positions": np.array(pos), "symbols": sym, "masses": np.array(mass),
            "fc": c.T @ sc["fc"] @ c / ncell, "dim": m}


def born_on_unit_cell(cell, run):
    """Z* of each unit-cell atom from the dielectric run (unit cell or any supercell of it):
    the average over the dielectric-run atoms that sit on that atom."""
    frac = run["positions"] @ run["lattice"] @ np.linalg.inv(cell["lattice"])
    n = len(cell["symbols"])
    born = np.zeros((n, 3, 3))
    count = np.zeros(n)
    for p, s, z in zip(frac, run["symbols"], run["born"]):
        k = find(p, s, cell["lattice"], cell["positions"], cell["symbols"])
        if k is None:
            raise ValueError(f"the atoms of {run['path']} do not sit on those of the phonon cell: "
                             "use the same structure (and origin) for both runs, or give --born")
        born[k] += z
        count[k] += 1
    if (count == 0).any():
        raise ValueError(f"{run['path']} has no Born charge for some atoms of the phonon cell")
    return born / count[:, None, None]
