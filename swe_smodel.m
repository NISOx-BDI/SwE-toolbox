function swe_smodel
% Launches batch window containing SwE batch module.
% =========================================================================
% This function prepares (if needed) and launches the batch system with a 
% job containing the batch module for the specification of
% data and design.
% =========================================================================
% FORMAT swe_smodel
% =========================================================================
% Written by Bryan Guillaume
% Version Info:  $Format:%ci$ $Format:%h$

    % Initiate a job 
    if isempty(spm_figure('FindWin','Graphics'))
        % SPM not running
        swe_jobman('initcfg')
    end
    
    % Launch the batch system with the SwE tab
    swe_batch
    % Add the specification module to it
    swe_jobman('interactive','','swe.smodel');

    return
    
end