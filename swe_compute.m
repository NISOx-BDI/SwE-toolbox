function swe_compute
% Launches batch window containing SwE batch module for running a design.
% =========================================================================
% This function prepares (if needed) and launches the batch system with a 
% job containing the batch module for the computation of a prespecified
% design.
% =========================================================================
% FORMAT swe_compute
% =========================================================================
% Written by Bryan Guillaume

    % Initiate a job 
    if isempty(spm_figure('FindWin','Graphics'))
        % SPM not running
        spm_jobman('initcfg')
    end
    
    % Launch the batch system with the SwE tab
    swe_batch
    % Add the specification module to it
    spm_jobman('interactive','','swe.compute');

    return
    
end