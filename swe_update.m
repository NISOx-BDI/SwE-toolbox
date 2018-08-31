% =========================================================================
% This function can be called with no input arguments to check whether or
% not the installed version of swe is up to date. It does this by looking
% on GitHub to find the latest version release number and comparing it to
% the version number stored in the local version of swe.
%
% Authors: Thomas Nichols, Tom Maullin (23/07/2018).
% =========================================================================

function swe_update()

    % Obtain swe version number.
    vswe=swe('ver');
    
    % Look to github for a newer version.
    url='https://github.com/NISOx-BDI/SwE-toolbox/releases/latest';
    [s,stat]=urlread(url);
    
    % Error if we couldn't contact github.
    if ~stat 
      error('Can''t contact GitHub');
    end
    
    % Look for latest swe version number.
    [tok,x]=regexp(s,'Version [0-9.]*','match','tokens');
    tok=tok{1};
    tok=strrep(tok,'Version ','');
    
    % Tell the user whether their version is the newest.
    if strcmp(tok,vswe)
      msg = 'Your version of swe is current.';
    else
      msg = ['A new version (' tok ') is available; please visit ' url ' to update your installation.'];
    end

    if ~nargout, fprintf([blanks(9) msg '\n']); end

end