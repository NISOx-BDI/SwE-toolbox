% ===============================================================
%
% This function was added for testing purposes and prevents the 
% spm_progress_bar displaying in octave. The progress bar can be 
% used in octave but it creates ascii art which spams the test
% logs. 
%
% Author: Tom Maullin (08/06/2018)
%
% ===============================================================
function swe_progress_bar(varargin)

	if ~exist('OCTAVE_VERSION','builtin')
       spm_progress_bar(varargin{:});
    end

end