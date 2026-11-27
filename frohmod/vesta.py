"""Draw a mode from an anime .ascii file (frohmod or phonopy) as arrows in VESTA.

    1. Open the POSCAR of the same calculation in VESTA and save it: File -> Save File -> structure.vesta
    2. frohmod ... --anime 6 --vesta structure.vesta        (or, for any anime file:
       python -m frohmod.vesta anime-6_q1_0_0.ascii structure.vesta anime-6.vesta)
    3. Open the new .vesta file in VESTA (arrow size and colour: Edit -> Vectors)

Use the POSCAR (not a CIF), so that VESTA keeps the atoms in the same order as the .ascii file.
"""
import sys

import numpy as np

longest = 1.0  # length of the longest arrow, in A


def ascii_to_vesta(ascii_file, vesta_in, vesta_out):
    """Write vesta_out: the structure of vesta_in with the first mode of ascii_file as black arrows."""
    with open(ascii_file) as f:
        lines = f.read().splitlines()

    # cell in V_Sim's frame (a along x, b in the xy plane)
    dxx, dyx, dyy = map(float, lines[1].split())
    dzx, dzy, dzz = map(float, lines[2].split())
    lat = np.array([[dxx, 0, 0], [dyx, dyy, 0], [dzx, dzy, dzz]])

    # displacements of the first mode in the file: the "#; ux; uy; uz; ..." lines
    u = []
    for l in lines:
        if l.startswith("#metaData") and u:
            break
        if l.startswith("#;"):
            u.append([float(x) for x in l[2:].rstrip(" \\").split(";")[:3]])
    u = np.array(u)

    # VESTA wants each arrow as components along a, b and c (in A)
    size = np.linalg.norm(u, axis=1)
    comp = (u @ np.linalg.inv(lat)) * np.linalg.norm(lat, axis=1) * longest / size.max()

    vectr, vectt = ["VECTR"], ["VECTT"]
    k = 0  # arrow number: 1, 2, 3, ... with no gaps, even when atoms are skipped
    for i, c in enumerate(comp, 1):
        if size[i - 1] < 0.01 * size.max():  # atom that hardly moves: no arrow
            continue
        k += 1
        vectr += [f"{k:4d} {c[0]:10.5f} {c[1]:10.5f} {c[2]:10.5f} 0", f"{i:5d}  0  0  0  0", " 0 0 0 0 0"]
        vectt.append(f"{k:4d}  0.200   0   0   0 1")  # thickness, then colour R G B: black
    vectr.append(" 0 0 0 0 0")
    vectt.append(" 0 0 0 0 0")

    # replace the (empty) VECTR and VECTT blocks of the saved file, which end where SPLAN starts
    with open(vesta_in) as f:
        v = f.read().splitlines()
    if "VECTR" not in v or "SPLAN" not in v:
        raise ValueError(f"{vesta_in} is not a .vesta file saved by VESTA")
    a, b = v.index("VECTR"), v.index("SPLAN")
    with open(vesta_out, "w") as f:
        f.write("\n".join(v[:a] + vectr + vectt + v[b:]) + "\n")
    return k


if __name__ == "__main__":
    n = ascii_to_vesta(*sys.argv[1:4])
    print(f"{sys.argv[3]}: {n} arrows")
