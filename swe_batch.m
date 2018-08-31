function swe_batch
% This function prepares and launches the batch system.
%
% This builds the whole tree for the various tools and their GUI at 
% the first call to this script.
% =========================================================================
% FORMAT bch = swe_batch
% =========================================================================
% Written by Bryan Guillaume

	% Check if SwE config tree is there
	if ~isstruct(cfg_util('tag2mod_cfg_id','swe.design'))
	    % SwE config tree
	    swe_gui = swe_cfg_batch;
	    % Adding SwE config tree to the SPM tools
	    cfg_util('addapp', swe_gui)
	end

	% Launching the batch system
	cfg_ui

end
