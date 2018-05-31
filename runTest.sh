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
octave --no-window-system --eval "moxunit_runtests()";