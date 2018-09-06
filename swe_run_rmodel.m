function swe_run_rmodel(varargin)
% Harvests job structure to run model stored in SwE.mat.
% =========================================================================
% FORMAT: swe_run_rmodel(job)
% -------------------------------------------------------------------------
% Inputs:
%  - job: harvested job data structure (see matlabbatch help)
% =========================================================================
% Written by Tom Maullin (05/09/2018)

% Job variable
% -------------------------------------------------------------------------
job   = varargin{1};

% Obtain SwE structure
% -------------------------------------------------------------------------
load(job.des{1});

% Compute model.
% -------------------------------------------------------------------------
cd(job.dir{1});
swe_cp(SwE);

end