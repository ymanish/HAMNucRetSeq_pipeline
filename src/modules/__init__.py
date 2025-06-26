import sys
from pathlib import Path

backend_dir = Path("~/pol/Projects/Codebase/Spermatogensis/backend").expanduser()
# backend_dir = Path("~/Spermatogensis_backend/backend").expanduser() ### Use it for the slrum cluster


# Verify base directory exists
if not backend_dir.exists():
    raise FileNotFoundError(f"Base directory not found: {backend_dir}")

# Add to Python path
nucfree_energy = backend_dir / "NucFreeEnergy"
sys.path.insert(0, str(nucfree_energy))


NUC_STATE_PATH = nucfree_energy / "methods" / "State" / "Nucleosome.state"
if not NUC_STATE_PATH.exists():
    raise FileNotFoundError(f"State file missing: {NUC_STATE_PATH}")

K_POSRESC_PATH = nucfree_energy / "MDParams" / "nuc_K_posresc_sym.npy"
if not K_POSRESC_PATH.exists():
    raise FileNotFoundError(f"Parameter file missing: {K_POSRESC_PATH}")

__all__ = ['NUC_STATE_PATH', 'K_POSRESC_PATH']

#### __all__ is a list of strings that explicitly declares which names (variables, functions, classes) should be considered "public" 
# when a user imports the module with from module import *.
#  If __all__ is defined, from module import * will only import the names listed in __all__.