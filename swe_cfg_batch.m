function swe = swe_cfg_batch
% This function builds the whole tree for the various tools and their GUI.
% =========================================================================
% FORMAT swe = swe_cfg_batch
% =========================================================================
% Written by Bryan Guillaume
    
    toolboxDir = spm_str_manip(mfilename('fullpath'), 'h');
    addpath(toolboxDir);
    addpath(fullfile(toolboxDir, 'test'));

	swe         = cfg_choice;
	swe.tag     = 'swe';
	swe.name    = 'SwE';
	swe.help    = {[...
	    'This is the batch interface for SwE, i.e. Sandwich Estimator '...
	    'Method for Neuroimaging Longitudinal Data Analysis.']
	                  }';
	swe.values  = {swe_config_smodel swe_config_rmodel swe_config_results};
	              
	%----------------------------------------------------------------------

end
