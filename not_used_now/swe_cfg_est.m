function swe_est = swe_cfg_est
% SwE Configuration file for Model Estimation


% $Id: spm_cfg_fmri_est.m 3753 2010-03-05 13:06:47Z guillaume $

% ---------------------------------------------------------------------
% spmmat Select SwE.mat
% ---------------------------------------------------------------------
swemat         = cfg_files;
swemat.tag     = 'swemat';
swemat.name    = 'Select SwE.mat';
swemat.help    = {
                  'Select the SwE.mat file that contains the design specification. '
                  'The directory containing this file is known as the input directory.'
}';
swemat.filter  = 'mat';
swemat.ufilter = '^SwE\.mat$';
swemat.num     = [1 1];

% ---------------------------------------------------------------------
% fmri_est Model estimation
% ---------------------------------------------------------------------
swe_est          = cfg_exbranch;
swe_est.tag      = 'swe_est';
swe_est.name     = 'Model estimation';
swe_est.val      = {swemat method };
swe_est.prog     = @swe_run_est;
swe_est.vout     = @vout_stats;

%==========================================================================
function dep = vout_stats(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'SwE.mat File';
dep(1).src_output = substruct('.','spmmat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
%dep(2)            = cfg_dep;
%dep(2).sname      = 'SPM Variable';
%dep(2).src_output = substruct('.','spmvar');
%dep(2).tgt_spec   = cfg_findspec({{'strtype','e'}});
if isfield(job.method, 'Classical')
    dep(2)            = cfg_dep;
    dep(2).sname      = 'Beta Images';
    dep(2).src_output = substruct('.','beta');
    dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    dep(3)            = cfg_dep;
    dep(3).sname      = 'Analysis Mask';
    dep(3).src_output = substruct('.','mask');
    dep(3).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    dep(4)            = cfg_dep;
    dep(4).sname      = 'ResMS Image';
    dep(4).src_output = substruct('.','resms');
    dep(4).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    % can't check whether auto-generated contrasts are generated this is
    % specified in input SwE.mat, not this job
end;
