function swe_batch
% Sandwich Estimator for Neuroimaging Longitudinal Data Analysis, SwE.
%
% This function prepares and launches the batch system.
% This builds the whole tree for the various tools and their GUI at the
% first call to this script.

% Written by Bryan Guillaume

persistent batch_initialize

if isempty(batch_initialize) || ~batch_initialize
    % SwE config tree
    swe_gui = swe_cfg_batch;
    % Adding SwE config tree to the SPM tools
    cfg_util('addapp', swe_gui)
    % No need to do it again for this session
    batch_initialize = 1;
end

% Launching the batch system
cfg_ui

return
