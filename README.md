# frohmod

**Which longitudinal optical (LO) phonon modes govern Fröhlich scattering, and which atoms contribute to their vibrations?**

<p align="center">
  <img src="assets/logo.svg" alt="frohmod" width="420">
</p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24598595.svg)](https://doi.org/10.5281/zenodo.24598595)

`frohmod` is a Python-based tool for performing mode-resolved analysis of Fröhlich (polar-optical phonon) scattering using phonon and dielectric properties obtained from VASP calculations.

The code incorporates the direction-dependent non-analytical correction to account for LO–TO splitting and decomposes the Fröhlich scattering into contributions from individual LO phonon modes and constituent elements. It identifies which LO modes dominate carrier–phonon scattering and quantifies the contribution of each atomic species to the corresponding phonon vibrations.

Unlike approaches that describe polar-optical scattering using a single effective phonon frequency, `frohmod` resolves the contributions of individual LO modes, providing a microscopic picture of the underlying scattering mechanism.

The code also generates V_Sim-compatible animations to visualize the atomic displacements associated with selected LO modes.

## Installation

Install directly from GitHub:

```bash
pip install git+https://github.com/shwetaphd/frohmod
```

This installs the required dependencies, NumPy and spglib.

For generating figures, matplotlib is also required:

```bash
pip install "frohmod[plot] @ git+https://github.com/shwetaphd/frohmod"
```

Alternatively:

```bash
pip install matplotlib
```

## Usage

The basic command (phonons along x, y and z) is:

```bash
frohmod --phonon vasprun.xml --plot
```

For a supercell calculation with a separate Born effective charge and dielectric calculation, along z only:

```bash
frohmod --phonon sc/vasprun.xml --dielectric-born unit/vasprun.xml --q 0 0 1 --anime 10 18
```

### Command-line options

| Option | Description |
|---|---|
| `--phonon` | `vasprun.xml` containing the phonon calculation (`IBRION = 6` with `NFREE = 2`, or `IBRION = 8`). Can be from a primitive cell or supercell calculation. |
| `--dim` | Supercell dimensions of the `--phonon` calculation, specified using phonopy's DIM format (e.g., `2 2 1`, or nine numbers). Default: found automatically from the structure. |
| `--dielectric-born` | `vasprun.xml` containing Born effective charges and the high-frequency dielectric tensor, calculated using `LCALCEPS` or `LEPSILON`. Can correspond to the unit cell or the same supercell. Default: the `--phonon` file. |
| `--q` | Cartesian direction of phonon propagation. Give it several times for several directions (e.g., `--q 1 0 0 --q 1 1 0`). Default: x, y and z. |
| `-T` | Temperature in K. Default: 300 K. |
| `--born` | Use an existing BORN file in phonopy format instead of generating one. |
| `--anime` | Mode numbers to animate, as printed in the `mode` column of the table. |
| `--vesta` | With `--anime`: also draw each mode as arrows in VESTA. Give the POSCAR of the same calculation, opened in VESTA and saved as a `.vesta` file. |
| `--plot` | Save the figures: `--plot modes` (modes.png, \|g\|² × n of every LO mode, split by element), `--plot elements` (elements.png, each element's share), or `--plot` alone for both. |


### Output files

The code generates the following files in the current working directory:

- `BORN`: Born effective charges and high-frequency dielectric tensor, unless an existing BORN file is supplied using `--born`.
- `modes.dat`: every mode of every direction: \|g\|², n, \|g\|² × n and each element's vibrational weight.
- `anime-<mode>_q<q>.ascii`: V_Sim-compatible animation files for the selected modes, one per direction.
- `anime-<mode>_q<q>.vesta`: the same modes as arrows in VESTA, with `--vesta`.
- `modes.png`, `elements.png`: the figures, with `--plot`.

## Output and interpretation

BiOBr (primitive cell, 300 K), phonons along x, y and z:

```text
q = (1.00, 0.00, 0.00)   T = 300 K   eps_inf(q) = 6.03   eps_0(q) = 34.78
cross-check: eps_0(q) from VASP's epsilon_ion = 34.90

  LO (THz) <- Gamma (THz)  mode  |g|^2 (meV)   n_LO  |g|^2 x n (meV)  Bi (meV)  Br (meV)   O (meV)
      3.74 <-        2.11     6         0.43   1.22             0.53      0.12      0.39      0.02
     11.44 <-        8.45    13         1.93   0.19             0.37      0.03      0.00      0.33
                              sum over all LO modes             0.89      0.15      0.39      0.36

(y: the same as x)

q = (0.00, 0.00, 1.00)   T = 300 K   eps_inf(q) = 5.72   eps_0(q) = 10.92
cross-check: eps_0(q) from VASP's epsilon_ion = 10.91

  LO (THz) <- Gamma (THz)  mode  |g|^2 (meV)   n_LO  |g|^2 x n (meV)  Bi (meV)  Br (meV)   O (meV)
      3.89 <-        3.13    10         0.39   1.16             0.46      0.11      0.34      0.01
     15.24 <-       13.71    18         1.08   0.10             0.10      0.01      0.00      0.10
                              sum over all LO modes             0.56      0.12      0.34      0.10

Frohlich coupling at 300 K from all LO modes along (1 0 0), (0 1 0), (0 0 1):  Bi 17.9 %   Br 47.4 %   O 34.7 %
```

### Understanding the output

| Column | Description |
|---|---|
| LO (THz) | Frequency of the LO mode for phonons propagating along q. |
| Gamma (THz) | The same vibration at Γ without the long-range field (no non-analytical correction). |
| mode | Its mode number at Γ, following phonopy's numbering (use it with `--anime`). |
| \|g\|² (meV) | Fröhlich coupling of the mode per phonon, (ħω/2)Δ. Independent of temperature. |
| n_LO | Number of phonons in the mode at temperature T (Bose–Einstein). |
| \|g\|² × n (meV) | Coupling at temperature T. |
| Bi, Br, O (meV) | Each element's part of \|g\|² × n: \|g\|² × n times the element's share of the mode's vibration. |
| Last line | Each element's share (%): its parts added over all LO modes of all directions, divided by the total. |

Modes below 3 % of a direction's total \|g\|² × n are not printed; all modes are kept in `modes.dat`. The last line depends on the directions given: along x only, BiOBr gives Bi 16.8 %, Br 43.3 %, O 39.9 %.

### Identifying and visualizing the dominant LO modes

The output table identifies the LO modes that contribute most strongly to Fröhlich scattering.

To visualize their atomic displacements:

1. Run `frohmod` to obtain the mode-resolved table.
2. Identify the LO modes with the largest \|g\|² × n.
3. Note their numbers in the `mode` column.
4. Run the calculation again using `--anime` with these mode numbers.

For example:

```bash
frohmod --phonon vasprun.xml --q 1 0 0 --anime 6 13
```

The generated `anime-<mode>_q<q>.ascii` files follow the same V_Sim format as phonopy's `anime.ascii`.

Importantly, these animations show each mode as it becomes along q, i.e. its LO eigenvector along the specified q direction. The LO displacement is given by:

```math
u = \frac{e}{\sqrt{M}}
```

Thus, unlike a conventional phonopy animation at Γ, these files visualize the atomic displacements associated with the direction-dependent LO modes responsible for Fröhlich scattering.

To draw a mode as arrows in VESTA (for a figure), open the POSCAR of the same calculation in VESTA, save it as `structure.vesta`, and add `--vesta`:

```bash
frohmod --phonon vasprun.xml --q 1 0 0 --anime 6 13 --vesta structure.vesta
```

This writes `anime-6_q1_0_0.vesta` and `anime-13_q1_0_0.vesta` next to the `.ascii` files. Any anime file (also phonopy's `anime.ascii`) can be converted on its own with:

```bash
python -m frohmod.vesta anime.ascii structure.vesta anime.vesta
```

The arrows are black and the longest is 2 Å; the arrow size and colour can be changed in VESTA under **Edit → Vectors**. Use the POSCAR rather than a CIF, so that VESTA keeps the atoms in the same order as the `.ascii` file.

## Figures and data files

`--plot` saves two figures for the run (needs matplotlib):

- `modes.png` (`--plot modes`): one bar per LO mode of every direction, height \|g\|² × n in meV, split by element.
- `elements.png` (`--plot elements`): each element's share (%) over all LO modes.

To compare several materials in one figure, run `frohmod` in each material's folder and then:

```bash
python -m frohmod.plot BiOBr/modes.dat Bi4O5Br2/modes.dat Bi3O4Br/modes.dat Bi5O7Br/modes.dat
```

This makes `modes.png` (one panel per material) and `elements.png` (one bar per material), labelled with the folder names.

`modes.dat` holds every mode of every direction, one block per direction:

| Column | Description |
|---|---|
| `qx`, `qy`, `qz` | Direction of the block. |
| `LO_THz`, `Gamma_THz`, `mode` | LO frequency, the frequency without the long-range field, and the mode number. |
| `g2_meV`, `n_LO`, `g2n_meV` | \|g\|², phonon number and \|g\|² × n. |
| `motion_<X>` | Fraction of the mode's vibrational energy carried by element X (0–1). |
| `part_<X>_meV` | \|g\|² × n × `motion_<X>`: element X's part of the mode. |

Element X's share is the sum of `part_<X>_meV` divided by the sum of `g2n_meV`.

## Methodology

The analysis is performed at the Brillouin-zone centre ($\Gamma$), for a specified phonon propagation direction $\hat{\mathbf q}$.

### 1. Modes at Γ without the long-range field

The phonon modes at Γ, without the macroscopic electric field, are obtained by diagonalizing the mass-weighted dynamical matrix:

```math
D_{k\alpha,k'\beta}
=
\frac{\Phi_{k\alpha,k'\beta}}
{\sqrt{M_k M_{k'}}}
```

where $\Phi$ is the force-constant matrix and $M_k$ is the atomic mass.

For supercell calculations, the force constants are first folded onto the unit cell.

### 2. Longitudinal optical (LO) modes

The LO modes are obtained by incorporating the long-range electrostatic contribution through the non-analytical correction to the dynamical matrix:

```math
\begin{aligned}
D^{\mathrm{NA}}_{k\alpha,k'\beta}
&=
\frac{4\pi e^2}{\Omega}
\\
&\quad \times
\frac{
(\hat{\mathbf q}\cdot Z^*_k)_\alpha
(\hat{\mathbf q}\cdot Z^*_{k'})_\beta
}{
(\hat{\mathbf q}\cdot\boldsymbol{\varepsilon}_\infty\cdot\hat{\mathbf q})
\sqrt{M_k M_{k'}}
}
\end{aligned}
```

where:

- $Z^*_k$ is the Born effective charge tensor of atom $k$.
- $\boldsymbol{\varepsilon}_\infty$ is the high-frequency dielectric tensor.
- $\Omega$ is the unit-cell volume.
- $\hat{\mathbf q}$ is the unit vector along the phonon propagation direction.

The corrected dynamical matrix is diagonalized to obtain the direction-dependent LO frequencies and eigenvectors.

### 3. Mode effective dipole

The effective mode dipole projected along the phonon propagation direction is calculated as:

```math
\hat{\mathbf q}\cdot\mathbf p_\nu
=
\hat{\mathbf q}\cdot
\sum_k
\frac{Z^*_k\mathbf e_{k\nu}}{\sqrt{M_k}}
```

where $\mathbf e_{k\nu}$ is the eigenvector of mode $\nu$.

This quantity describes the coupling of the phonon mode to the macroscopic electric field.

### 4. Fröhlich coupling per phonon

Each LO mode's share of the polar screening:

```math
\Delta_\nu
=
\frac{4\pi e^2}{\Omega}
\frac{(\hat{\mathbf q}\cdot\mathbf p_\nu)^2}
{\omega_\nu^2\,(\hat{\mathbf q}\cdot\boldsymbol{\varepsilon}_\infty\cdot\hat{\mathbf q})^2},
\qquad
\sum_\nu \Delta_\nu
=
\frac{1}{\varepsilon_\infty(\hat{\mathbf q})}-\frac{1}{\varepsilon_0(\hat{\mathbf q})}
```

The coupling of the mode is AMSET's polar-optical matrix element (Ganose et al., 2021), one LO mode at a time, in the $q\rightarrow 0$ limit of Verdi and Giustino (2015):

```math
|g_\nu|^2
=
\frac{\hbar\omega_\nu}{2}\,\Delta_\nu
```

where $\omega_\nu$ is the frequency of the LO mode. The electronic overlap and the factor $4\pi e^2/(\Omega q^2)$ are the same for every mode and are left out.

### 5. Phonon occupation

The thermal occupation of each LO phonon mode at temperature $T$ is calculated using the Bose–Einstein distribution:

```math
n_\nu
=
\frac{1}{
\exp\left(\frac{\hbar\omega_\nu}{k_{\mathrm B}T}\right)-1
}
```

### 6. Coupling at temperature T

```math
G_\nu = |g_\nu|^2\, n_\nu \quad \text{(meV)}
```

### 7. Element-resolved contributions

The fraction of the vibrational energy of mode $\nu$ carried by element $X$:

```math
P_\nu(X)
=
\sum_{k\in X}
|\mathbf e_{k\nu}|^2
```

The share of element $X$, over all LO modes $\nu$ of all directions $\hat{\mathbf q}$:

```math
C_X
=
100\,
\frac{\sum_{\hat{\mathbf q}}\sum_\nu G_\nu P_\nu(X)}
{\sum_{\hat{\mathbf q}}\sum_\nu G_\nu}
```

## Notes and limitations

- \|g\|² × n is **not an absolute scattering rate**: only the long-range Fröhlich coupling in the $q\rightarrow 0$ limit is considered, without carrier mass, band structure or wavefunction overlap. It is the same for electrons and holes.
- The element shares depend on the directions given with `--q`; state them when reporting the shares.
- The `--phonon` and `--dielectric-born` calculations must use the same structure and origin. A rigid translation between the two structures is not automatically detected. In such cases, provide an appropriate BORN file using `--born`.

## Code structure

| File | Description |
|---|---|
| `frohmod/read.py` | Reads and extracts phonon, dielectric, and Born effective charge information from VASP `vasprun.xml` files. |
| `frohmod/supercell.py` | Finds the supercell size and maps supercell force constants and Born effective charges onto the unit cell. |
| `frohmod/born.py` | Reads and writes phonopy-compatible BORN files using crystal symmetry. |
| `frohmod/analysis.py` | Implements the mode-resolved Fröhlich analysis (steps 1–7) and the element shares. |
| `frohmod/output.py` | Generates output tables, `.dat` files, figures, and V_Sim-compatible animations. |
| `frohmod/cli.py` | Implements the `frohmod` command-line interface. |
| `frohmod/plot.py` | Draws the two figures from `modes.dat` files (`--plot`, or several materials with `python -m frohmod.plot`). |
| `frohmod/vesta.py` | Draws an `anime` mode as arrows in VESTA (`--vesta`). |

## Author

**Shweta Choudhary**

Department of Chemistry  
Indian Institute of Technology Roorkee, India 

Email: shwetachoudharyphd@gmail.com

ORCID: [0000-0002-6449-9459](https://orcid.org/0000-0002-6449-9459)

## Acknowledgements

Shweta Choudhary acknowledges the financial support provided by the Prime Minister's Research Fellowship (PMRF), Government of India.

Claude (Anthropic) was used in developing the code.

## Citing

If you use `frohmod` in your research, please cite this repository using the **Cite this repository** option on GitHub, based on the information provided in `CITATION.cff`.

Please also cite the relevant methodological and computational references:

1. C. Verdi and F. Giustino, *Physical Review Letters* **115**, 176401 (2015).
2. X. Gonze and C. Lee, *Physical Review B* **55**, 10355 (1997).
3. A. Togo and I. Tanaka, *Scripta Materialia* **108**, 1 (2015).
4. A. M. Ganose et al., *Nature Communications* **12**, 2222 (2021).

## License

This project is distributed under the MIT License.
