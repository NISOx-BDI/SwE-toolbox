function results = swe_config_results
% This displays the results from an SwE model.
% =========================================================================
% FORMAT results = swe_config_results
% =========================================================================
% Written by Tom Maullin (05/09/2018)
% Version Info:  $Format:%ci$ $Format:%h$

% -------------------------------------------------------------------------
% Results Directory
% -------------------------------------------------------------------------
dis         = cfg_files;
dis.tag     = 'dis';
dis.name    = 'SwE File';
dis.help    = {' '
    'Select an ''SwE.mat'' file containing the computed model to be displayed.'};
dis.filter = 'mat';
dis.ufilter = '.*';
dis.num     = [1 1];

% -------------------------------------------------------------------------
% compute design
% -------------------------------------------------------------------------
results        = cfg_exbranch;
results.tag    = 'results';
results.name   = 'Results Display';
results.val    = {dis};
results.help   = {' '
                 'Module of the SwE toolbox allowing the display of a previously computed model.'};
results.prog   = @swe_run_results;