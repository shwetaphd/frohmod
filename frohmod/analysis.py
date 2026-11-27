"""The mode-resolved Frohlich analysis at Gamma"""
import numpy as np

thz = 15.633302      # sqrt(eV / A^2 / amu) -> THz
coulomb = 14.399645  # e^2 / (4 pi eps0) in eV A
kb = 0.08617333      # Boltzmann constant, meV / K
h = 4.135667696      # Planck constant, meV / THz


def freq(w2):
    """Frequencies in THz from eigenvalues w^2; imaginary modes come out negative."""
    return np.sign(w2) * np.sqrt(np.abs(w2)) * thz


def analyse(cell, born, eps, direction, temperature=300.0):
    """cell: unit cell with force constants (dict from read.py / supercell.py),
    born: Z* of every atom (n x 3 x 3), eps: eps_inf (3 x 3), direction: q (Cartesian).
    Returns a dict with the LO frequencies, the frequencies without the long-range field, the coupling |g|^2, phonon number n and |g|^2 x n (meV)
    of every LO mode, and each element's share of every mode."""
    q = np.array(direction, float)
    q /= np.linalg.norm(q)
    s, masses = cell["symbols"], cell["masses"]
    n = len(s)
    m = np.repeat(masses, 3)
    vol = abs(np.linalg.det(cell["lattice"]))
    eq = q @ eps @ q  # eps_inf along q

    # 1. modes at Gamma without the long-range field (no NAC)
    dyn = cell["fc"] / np.sqrt(np.outer(m, m))
    w2_gamma, e_gamma = np.linalg.eigh(dyn)

    # 2. LO modes: add the field made by the charges moving along q (non-analytical correction)
    born_sum = np.abs(born.sum(axis=0)).max()  # should be ~0
    qz = np.array([q @ z for z in born - born.mean(axis=0)])
    v = (qz / np.sqrt(masses)[:, None]).ravel()  # (q . Z*_k) / sqrt(M_k) for every atom k
    w2_lo, e_lo = np.linalg.eigh(dyn + 4 * np.pi * coulomb / vol * np.outer(v, v) / eq)
    f_gamma, f_lo = freq(w2_gamma), freq(w2_lo)

    # 3. dipole of each mode along q: q . sum_k Z*_k e_k / sqrt(M_k), which is just v @ e
    # 4. Delta: the mode's share of 1/eps_inf - 1/eps_0 (all Delta add up to it), and its
    #    Frohlich coupling |g|^2 = (hbar w / 2) Delta in meV (AMSET's matrix element, one LO mode at a time;
    #    the electron overlap and 4 pi e^2 / (Omega q^2) are the same for every mode and left out)
    # 5. phonons present at T (Bose-Einstein), and the coupling at T: |g|^2 x n
    ok = f_lo > 0.1  # leave out the acoustic modes
    delta = np.zeros(3 * n)
    nph = np.zeros(3 * n)
    delta[ok] = 4 * np.pi * coulomb / vol * (v @ e_lo[:, ok])**2 / (w2_lo[ok] * eq**2)
    nph[ok] = 1 / np.expm1(h * f_lo[ok] / (kb * temperature))
    g2 = 0.5 * h * np.clip(f_lo, 0, None) * delta
    g2n = g2 * nph

    # 6. who moves: each element's share of a mode's vibrational energy (sum of |e|^2 over its atoms),
    #    and its part of the mode's |g|^2 x n (meV)
    a = (e_lo**2).reshape(n, 3, -1).sum(axis=1)
    elements = list(dict.fromkeys(s))
    motion = {}
    part = {}
    for el in elements:
        motion[el] = a[[x == el for x in s]].sum(axis=0)
        part[el] = g2n * motion[el]

    # which no-NAC mode each LO mode comes from (largest overlap)
    overlap = (e_gamma.T @ e_lo)**2

    # 7. eps_0 along q rebuilt from the modes without NAC (compared with VASP's epsilon_ion in the table)
    ok = f_gamma > 0.1
    eps0 = eq + 4 * np.pi * coulomb / vol * ((v @ e_gamma[:, ok])**2 / w2_gamma[ok]).sum()

    return {"q": q, "temperature": temperature, "elements": elements, "f_lo": f_lo, "f_gamma": f_gamma, "e_lo": e_lo,
            "overlap": overlap, "gamma_mode": overlap.argmax(axis=0), "delta": delta, "g2": g2, "n": nph, "g2n": g2n,
            "motion": motion, "part": part, "eps_q": eq, "eps0_q": eps0, "born_sum": born_sum,
            "n_imaginary": int((f_gamma < -0.1).sum())}


def element_share(results):
    """Each element's share (%) of |g|^2 x n, over all LO modes of all directions together."""
    els = results[0]["elements"]
    total = {el: sum(r["part"][el].sum() for r in results) for el in els}
    return {el: 100 * total[el] / sum(total.values()) for el in els}
