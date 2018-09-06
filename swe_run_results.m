function swe_run_results(varargin)
% Harvests job structure to display results stored in SwE.mat.
% =========================================================================
% FORMAT: swe_run_results(job)
% -------------------------------------------------------------------------
% Inputs:
%  - job: harvested job data structure (see matlabbatch help)
% =========================================================================
% Written by Tom Maullin (05/09/2018)

% Job variable
% -------------------------------------------------------------------------
job   = varargin{1};

% Change directory to directory containing SwE structure
% -------------------------------------------------------------------------
cd(fileparts(job.dis{1}));

% Display results.
% -------------------------------------------------------------------------
swe_results_ui; %%ISSUE: DISPLAY OUTPUT NOT IN USER PREFERENCE FOLDER

end