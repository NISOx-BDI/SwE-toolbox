function swe_results
% Launches batch window containing SwE batch module.
% =========================================================================
% This function prepares (if needed) and launches the batch system with a 
% job containing the batch module for the specification of
% data and design.
% =========================================================================
% FORMAT swe_results
% =========================================================================
% Written by Tom Maullin (05/09/2018)
% Version Info:  $Format:%ci$ $Format:%h$

    % Initiate a job 
    if isempty(spm_figure('FindWin','Graphics'))
        % SPM not running
        spm_jobman('initcfg')
    end
    
    % Launch the batch system with the SwE tab
    swe_batch
    % Add the specification module to it
    spm_jobman('interactive','','swe.results');

    return
    
end