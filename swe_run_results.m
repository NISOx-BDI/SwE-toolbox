function swe_run_results(varargin)
% Harvests job structure to display results stored in SwE.mat.
% =========================================================================
% FORMAT: swe_run_results(job)
% -------------------------------------------------------------------------
% Inputs:
%  - job: harvested job data structure (see matlabbatch help)
% -------------------------------------------------------------------------
% Outputs:
%  - out: filename of saved data structure.
% =========================================================================
% Written by Tom Maullin (05/09/2018)

% Job variable
% -------------------------------------------------------------------------
job   = varargin{1};

% Obtain SwE structure
% -------------------------------------------------------------------------
load(job.res{1});

% Display results.
% -------------------------------------------------------------------------
swe_cp(SwE);

end