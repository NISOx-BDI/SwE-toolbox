function results = swe_cfg_results
% This displays the results from an SwE model.
% =========================================================================
% FORMAT modign = swe_cfg_results
% =========================================================================
% Written by Tom Maullin (05/09/2018)

% -------------------------------------------------------------------------
% dir Directory
% -------------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {' '
    'Select a directory where the output images will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% -------------------------------------------------------------------------
% dir Directory
% -------------------------------------------------------------------------
mod         = cfg_files;
mod.tag     = 'mod';
mod.name    = 'Model File';
mod.help    = {' '
    'Select an ''SwE.mat'' file containing the computed model to be displayed.'};
mod.filter = 'mat';
mod.ufilter = '.*';
mod.num     = [1 1];

% -------------------------------------------------------------------------
% compute modign
% -------------------------------------------------------------------------
results        = cfg_exbranch;
results.tag    = 'results';
results.name   = 'Results Display';
results.val    = {dir mod};
results.help   = {' '
                 'Module of the SwE toolbox allowing the display of a previously computed model.'};
results.prog   = @swe_run_results;%%%%%should be run results