function swe = swe_cfg_batch
% Sandwich Estimator Method for Neuroimaging Longitudinal Data Analysis, SwE.
% SwE configuration file
% This builds the whole tree for the various tools and their GUI.

% Written by Bryan Guillaume

swe         = cfg_choice;
swe.tag     = 'swe';
swe.name    = 'SwE';
swe.help    = {[...
    'This is the batch interface for SwE, i.e. Sandwich Estimator '...
    'Method for Neuroimaging Longitudinal Data Analysis.']
                  }';
swe.values  = {swe_cfg_design};
              
%------------------------------------------------------------------------

return