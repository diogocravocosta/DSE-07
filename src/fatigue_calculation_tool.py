import numpy as np

from pathlib import Path
import sys

current_dir = Path(__file__).parent
sys.path.append(str(current_dir.parent))

from data import constants as cs


print(cs.g_0)

def thermal_stress(delta_T, young_mod, thermal_expansion_coeff):
    return delta_T * young_mod * thermal_expansion_coeff



