function swe = swe_cfg_batch
% This function builds the whole tree for the various tools and their GUI.
% =========================================================================
% FORMAT swe = swe_cfg_batch
% =========================================================================
% Written by Bryan Guillaume

    addpath

	swe         = cfg_choice;
	swe.tag     = 'swe';
	swe.name    = 'SwE';
	swe.help    = {[...
	    'This is the batch interface for SwE, i.e. Sandwich Estimator '...
	    'Method for Neuroimaging Longitudinal Data Analysis.']
	                  }';
	swe.values  = {swe_leaf_smodel swe_leaf_rmodel swe_leaf_results};
	              
	%----------------------------------------------------------------------

end
