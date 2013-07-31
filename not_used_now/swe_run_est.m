function out = swe_run_fmri_est(job)
% Estimate parameters of a specified model
% SwE job execution function
% takes a harvested job data structure and call SPM ans SwE functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________


%-Load SwE.mat file
%--------------------------------------------------------------------------
SwE = [];
load(job.swemat{:});
out.swemat = job.swemat;

%-Move to the directory where the SwE.mat file is
%--------------------------------------------------------------------------
original_dir = pwd;
cd(fileparts(job.spmmat{:}));



%out.spmvar = SPM;
cd(original_dir);
fprintf('Done\n')
return