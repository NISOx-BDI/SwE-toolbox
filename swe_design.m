function swe_design
% Sandwich Estimator for Neuroimaging Longitudinal Data Analysis, SwE.
%
% This function prepares (if needed) and launches the batch system with a job containing the batch module for the specification of
% data and design.

% Written by Bryan Guillaume

% Initiate a job 
spm_jobman('initcfg')
% Launch the batch system with the SwE tab
swe_batch
% Add the specification module to it
spm_jobman('interactive','','swe.design');

return