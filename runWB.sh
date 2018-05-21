mv code/spm12/\@file_array/private/mat2file.mex .
rm code/spm12/\@file_array/private/mat2file.* 
mv ./mat2file.mex code/spm12/\@file_array/private/
cd swe
octave --no-window-system --eval "addpath(genpath(pwd));runTest();load('test/data/design.mat');swe_run_design(design);load('SwE.mat');swe_cp_WB(SwE);"

