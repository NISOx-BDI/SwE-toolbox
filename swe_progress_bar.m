function swe_progress_bar(varargin)
% Opens the progress bar, if not in octave.
% =========================================================================
% This function was added for testing purposes and prevents the
% spm_progress_bar displaying in octave. The progress bar can be used in
% octave but it creates ascii art which spams the test logs.
% =========================================================================
% FORMAT swe_progress_bar('Set',value)
% Set the height of the bar itself.
%
% FORMAT swe_progress_bar('Set','xlabel',xlabel)
% FORMAT swe_progress_bar('Set','ylabel',ylabel)
% Set the progress bar labels.
%
% FORMAT swe_progress_bar('Set','height',height)
% Set the height of the progress bar.
%
% FORMAT swe_progress_bar('Clear')
% Clear the 'Interactive' window.
% =========================================================================
% Author: Tom Maullin (08/06/2018)
% Version Info:  $Format:%ci$ $Format:%h$
% =========================================================================

    if ~exist('OCTAVE_VERSION','builtin')
        spm_progress_bar(varargin{:});
    end

end
