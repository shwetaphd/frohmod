"""frohmod: mode-resolved analysis of Frohlich scattering from VASP phonon calculations."""
from .analysis import analyse, coulomb, element_share, thz
from .cli import main
from .read import read_vasprun

__version__ = "0.1.0"
