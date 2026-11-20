#!/bin/bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Run E4 (local V-USRP + local CHEM) in WSL. Uses the 'vusrp' env (UHD with VUSRP, Python API).
export MAMBA_ROOT_PREFIX=~/mm/root
eval "$(~/mm/bin/micromamba shell hook -s bash)"
micromamba activate vusrp
PFX=~/work/vusrp/install
export LD_LIBRARY_PATH=$PFX/lib:$CONDA_PREFIX/lib
export PYTHONPATH=$(ls -d $PFX/lib/python3*/site-packages 2>/dev/null | head -1):$PYTHONPATH
cd "$(dirname "$0")"
python -u run_e4.py "$@" 2>&1
