"""Read vasprun.xml: the structure and masses, and whatever the run holds of force constants,
Born charges and dielectric tensors. Anything the run doesn't have comes back as None."""
import xml.etree.ElementTree as ET

import numpy as np

from .analysis import thz


def rows(x):
    """The numbers in the <v> rows of a vasprun.xml array, as a 2D array."""
    return np.array([[float(a) for a in v.text.split()] for v in x.findall("v")])


def read_vasprun(path):
    """Return a dict with path, lattice, positions (fractional), symbols, masses,
    fc (force constants, eV/A^2), born (Z* of every atom), eps (eps_inf) and eps_ion."""
    try:
        root = ET.parse(path).getroot()
    except FileNotFoundError:
        raise ValueError(f"{path} not found")
    except ET.ParseError:
        raise ValueError(f"{path} is incomplete: did the VASP run finish?")

    # atoms and their masses (POMASS of each atom type)
    atoms = root.find("atominfo/array[@name='atoms']/set")
    if atoms is None:
        raise ValueError(f"{path} is not a VASP vasprun.xml")
    symbols = [r[0].text.strip() for r in atoms]
    mass = [float(r[2].text) for r in root.find("atominfo/array[@name='atomtypes']/set")]
    masses = np.array([mass[int(r[1].text) - 1] for r in atoms])

    # the undisplaced structure (with IBRION = 6, the final structure and CONTCAR are the last displaced one)
    s = root.find("structure[@name='initialpos']")
    run = {"path": str(path), "symbols": symbols, "masses": masses,
           "lattice": rows(s.find("crystal/varray[@name='basis']")),
           "positions": rows(s.find("varray[@name='positions']")),
           "fc": None, "born": None, "eps": None, "eps_ion": None}

    # force constants: VASP stores -Phi/sqrt(M M'), in THz^2 (VASP 6) or in eV/A^2/amu (older versions)
    h = root.findall(".//dynmat/varray[@name='hessian']")
    if h:
        h = rows(h[-1])
        unit = root.find(".//dynmat/i[@name='unit']")
        if unit is not None and "THz" in unit.text:
            h = h / thz**2
        m = np.repeat(masses, 3)
        fc = -h * np.sqrt(np.outer(m, m))
        run["fc"] = (fc + fc.T) / 2  # make it exactly symmetric

    # Born charges, one 3x3 tensor per atom
    z = root.findall(".//array[@name='born_charges']")
    if z:
        run["born"] = np.array([rows(a) for a in z[-1].findall("set")])

    # eps_inf: LEPSILON writes 'epsilon', LCALCEPS writes 'epsilon_scf'
    e = root.findall(".//varray[@name='epsilon']") or root.findall(".//varray[@name='epsilon_scf']")
    if e:
        run["eps"] = rows(e[-1])

    # ionic part of eps_0, only used to cross-check the result
    e = root.findall(".//varray[@name='epsilon_ion']")
    if e:
        run["eps_ion"] = rows(e[-1])
    return run
