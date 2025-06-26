# src/config/path.py
# Author: MY
# Created on: 2025-04-02


from pathlib import Path
import datetime


######LOCAL DIR#######

SRC_DIR = Path(__file__).parent.parent.parent
OUT_DIR = Path(__file__).parent.parent.parent


######CLUSTER DIR#######

# SRC_DIR = Path('/home/pol_schiessel/maya620d/NucRetention_pipeline')
# OUT_DIR = Path('/group/pol_schiessel/Manish/NucRetention_pipeline')


print(f"Source Directory: {SRC_DIR}")
print(f"Output Directory: {OUT_DIR}")


PARAM_DIR = SRC_DIR / 'parameters'
TEST_DIR = SRC_DIR / 'tests'
EXAMPLES_DIR = SRC_DIR / 'examples'


DATA_DIR = OUT_DIR / 'data'
RESULTS_DIR = OUT_DIR / 'output'



