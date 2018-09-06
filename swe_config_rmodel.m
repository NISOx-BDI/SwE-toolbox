function rmodel = swe_config_rmodel
% This takes in a SwE.mat file containing a design and runs the model.
% =========================================================================
% FORMAT design = swe_config_rmodel
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
des         = cfg_files;
des.tag     = 'des';
des.name    = 'SwE File';
des.help    = {' '
    'Select an ''SwE.mat'' file containing the model to be run.'};
des.filter = 'mat';
des.ufilter = '.*';
des.num     = [1 1];

% ---------------------------------------------------------------------
% compute Design
% ---------------------------------------------------------------------
rmodel        = cfg_exbranch;
rmodel.tag    = 'rmodel';
rmodel.name   = 'Run Model';
rmodel.val    = {dir des};
rmodel.help   = {' '
                 'Module of the SwE toolbox allowing the computing of a previously specified model.'};
rmodel.prog   = @swe_run_rmodel;
