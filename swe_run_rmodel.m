function swe_run_rmodel(varargin)
% Harvests job structure to run model stored in SwE.mat.
% =========================================================================
% FORMAT: swe_run_rmodel(job)
% -------------------------------------------------------------------------
% Inputs:
%  - job: harvested job data structure (see matlabbatch help)
% =========================================================================
% Written by Tom Maullin (05/09/2018)
% Version Info:  $Format:%ci$ $Format:%h$

% Job variable
% -------------------------------------------------------------------------
job   = varargin{1};

% Obtain SwE structure
% -------------------------------------------------------------------------
load(job.des{1});

% Check compatability
% -------------------------------------------------------------------------
if isfield(SwE, 'ver')
    swe_checkCompat(SwE.ver, swe('ver'))
else
    warning(['SwE.mat file has no recorded version number. If this `.mat` ',...
             'file was created with SwE v1.2.11 or less, running this anal',...
             'ysis will likely run into errors. To avoid this, please re-e',...
             'nter the job specification in the batch window.']);
end

% Compute model.
% -------------------------------------------------------------------------
cd(SwE.swd);
swe_cp(SwE);

end
