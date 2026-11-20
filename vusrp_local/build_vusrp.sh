#!/bin/bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Build the public V-USRP (anilgurses/VUSRP, a UHD fork) in WSL, in its own micromamba env and prefix.
set -e
export MAMBA_ROOT_PREFIX=~/mm/root
eval "$(~/mm/bin/micromamba shell hook -s bash)"
if ! micromamba env list | grep -q " vusrp "; then
  micromamba create -y -n vusrp -c conda-forge python=3.11 cmake ninja cxx-compiler c-compiler pkg-config \
    libboost-devel yaml-cpp mako ruamel.yaml numpy pybind11 libusb requests 2>&1 | tail -3
fi
micromamba install -y -n vusrp -c conda-forge eigen "libboost-devel=1.84" 2>&1 | tail -2
micromamba activate vusrp
rm -rf ~/work/vusrp/build
SRC=~/work/vusrp/VUSRP/host; BLD=~/work/vusrp/build; PFX=~/work/vusrp/install
mkdir -p $BLD && cd $BLD
cmake $SRC -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PFX -DCMAKE_PREFIX_PATH=$CONDA_PREFIX "-DCMAKE_CXX_FLAGS=-include cstdint" \
  -DENABLE_VUSRP=ON -DENABLE_PYTHON_API=ON -DENABLE_EXAMPLES=ON -DENABLE_UTILS=ON -DENABLE_TESTS=OFF \
  -DENABLE_MANUAL=OFF -DENABLE_DOXYGEN=OFF -DENABLE_MAN_PAGES=OFF -DENABLE_USB=OFF -DENABLE_B100=OFF -DENABLE_B200=OFF \
  -DENABLE_X300=OFF -DENABLE_N300=OFF -DENABLE_N320=OFF -DENABLE_E300=OFF -DENABLE_E320=OFF -DENABLE_X400=OFF \
  -DENABLE_USRP1=OFF -DENABLE_USRP2=OFF -DENABLE_OCTOCLOCK=OFF -DENABLE_MPMD=ON -DENABLE_SIM=OFF > cmake.log 2>&1 || { grep -iE "error|vusrp|not found" cmake.log | head -30; exit 1; }
grep -iE "VUSRP|Building version" cmake.log | head -8
ninja -j8 > ninja.log 2>&1 || { grep -E "error|FAILED" ninja.log | head -30; exit 1; }
ninja install > install.log 2>&1
ls $PFX/bin | head; ls $PFX/lib/uhd/examples | head -40
echo VUSRP_BUILD_DONE
