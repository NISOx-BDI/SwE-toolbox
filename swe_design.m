function swe_design
% This function prepares (if needed) and launches the batch system with a 
% job containing the batch module for the specification of
% data and design.
% =========================================================================
% FORMAT swe_design
% =========================================================================
% Written by Bryan Guillaume

    % Initiate a job 
    spm_jobman('initcfg')
    % Launch the batch system with the SwE tab
    swe_batch
    % Add the specification module to it
    spm_jobman('interactive','','swe.design');

    return
    
end