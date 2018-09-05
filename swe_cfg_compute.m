function design = swe_cfg_compute
% This takes in the SwE.mat data and computes the design.
% =========================================================================
% FORMAT design = swe_cfg_compute
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
des.name    = 'Design File';
des.help    = {' '
    'Select an ''SwE.mat'' file containing the design to be compute.'};
des.filter = 'mat';
des.ufilter = '.*';
des.num     = [1 1];

% ---------------------------------------------------------------------
% compute Design
% ---------------------------------------------------------------------
design        = cfg_exbranch;
design.tag    = 'compute';
design.name   = 'Compute Design';
design.val    = {dir des};
design.help   = {' '
                 'Module of the SwE toolbox allowing the computing of a specified design.'};
design.prog   = @swe_cp;
design.vout   = @vout_data;
