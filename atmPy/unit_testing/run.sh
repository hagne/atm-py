#!/bin/bash

PYTHONPATH="/Users/htelg/prog/hagpack/:$PYTHONPATH"
export PYTHONPATH
PYTHONPATH="/Users/htelg/projecte/POPS/prog/skriptchen:$PYTHONPATH"
export PYTHONPATH
PYTHONPATH="/Users/htelg/projecte/POPS/mie_scattering/bhmie-f/:$PYTHONPATH"
export PYTHONPATH
PYTHONPATH="/Users/htelg/prog/atm-py/:$PYTHONPATH"
export PYTHONPATH
PYTHONPATH="/Users/htelg/prog/:$PYTHONPATH"
export PYTHONPATH
export PATH=/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin:$PATH
export PATH=/Users/htelg/prog/skriptchen/:$PATH
export PATH=/Users/htelg/bin:$PATH
export PATH=/opt/local/bin:$PATH

/Users/htelg/prog/atm-py/atmPy/unit_testing/run_nose_tests.py