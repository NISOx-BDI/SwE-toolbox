#!/bin/bash
# =======================================================================
# This file contains the setup for the SwE-toolbox CI travis tests which
# MUST be done inside the docker container.
#
# Author: Tom Maullin (31/05/2018) 
# =======================================================================
# Octave-SPM bug fix
mv code/spm12/\@file_array/private/mat2file.mex .
rm code/spm12/\@file_array/private/mat2file.* 
mv ./mat2file.mex code/spm12/\@file_array/private/
cd swe

# Testing install
ls
cd ./test/MOxUnit
make install-octave
cd ..

# Octave commands
testresult=$(octave --no-window-system --eval "moxunit_runtests();")

echo "$testresult"

# If the tests failed, we need to let Travis know.
if [[ $testresult = *"FAILED"* ]]; then
  exit 1
else
  exit 0
fi