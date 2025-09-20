# src/config/path.py
# Author: MY
# Created on: 2025-04-02


from pathlib import Path
import datetime


######LOCAL DIR#######

SRC_DIR = Path(__file__).parent.parent.parent
OUT_DIR = Path(__file__).parent.parent.parent


######CLUSTER DIR#######

# SRC_DIR = Path('/home/pol_schiessel/maya620d/HAMNucRetSeq_pipeline')
# OUT_DIR = Path('/group/pol_schiessel/Manish/HAMNucRetSeq_pipeline')


print(f"Source Directory: {SRC_DIR}")
print(f"Output Directory: {OUT_DIR}")


PARAM_DIR = SRC_DIR / 'parameters'
TEST_DIR = SRC_DIR / 'tests'
EXAMPLES_DATA = SRC_DIR / 'example_data'
EXAMPLES_OUT = SRC_DIR / 'example_output'


DATA_DIR = OUT_DIR / 'data'
RESULTS_DIR = OUT_DIR / 'output'
NUCMONTE_DIR = OUT_DIR / 'nucmontecarlo_data' ### The data was copied from the NucMonteCarlo repository



