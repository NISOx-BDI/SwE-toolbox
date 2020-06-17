function swe_checkCompat(matVer, tbVer)
% Checks for compatibility between SwE versions and errors if there is a
% compatibility issue.
% =========================================================================
% FORMAT: swe_checkCompat(matVer, tbVer)
% -------------------------------------------------------------------------
% Inputs: 
%   - matVer: version of SwE recorded in the `.mat` file.
%   - tbVer:  version of SwE toolbox currently being run.
% =========================================================================
% Author: Tom Maullin (09/11/2018)
% Version Info:  $Format:%ci$ $Format:%h$

    if isequal(matVer,tbVer)
       return; 
    end
    
    % The below hashmap records the earliest compatible version for each
    % release of the SwE toolbox. I.e. when making a new release, say you
    % are releasing version "y.y.y", please set earliestCompatVer("y.y.y")
    % equal to "x.x.x" where "x.x.x" is the oldest version of the toolbox
    % which "y.y.y" can accept `SwE.mat` files from.
    earliestCompatVer = containers.Map();
    
    % These versions did not record version numbers in the SwE.mat file and
    % therefore cannot be checked in the same way. However, none of these
    % should be compatible with version 2.0.1 (the version in which this
    % was released) or anything further.
    earliestCompatVer('1.0') = 'NaN.NaN.NaN';
    earliestCompatVer('1.1') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.1') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.2') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.3') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.4') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.5') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.6') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.7') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.8') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.9') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.10') = 'NaN.NaN.NaN';
    earliestCompatVer('1.2.11') = 'NaN.NaN.NaN';
    
    % Record earliest compatible versions.
    earliestCompatVer('2.0.0') = '2.0.0';
    earliestCompatVer('2.1.0') = '2.0.0';
    earliestCompatVer('2.1.1') = '2.0.0';
    earliestCompatVer('2.2.0.rc') = '2.0.0';
    earliestCompatVer('2.2.0.rc2') = '2.0.0';
    earliestCompatVer('2.2.0') = '2.0.0';
    earliestCompatVer('2.2.1') = '2.0.0';

    % The below line works out the latest compatible version from the
    % earliest compatible versions. This code is now redundant but may be
    % useful in future so has been left in place. Tom Maullin (09/11/2018)
    %latestCompatVer = latComVer(earliestCompatVer);
    
    % Check if the `.mat` version is compatible with this version.
    if ~strcmp(earliestCompatVer(matVer), earliestCompatVer(tbVer)) || ...
            strcmp(earliestCompatVer(matVer), 'NaN.NaN.NaN')
        error(['The SwE version used to create this `.mat` file (version ',...
               matVer, ') is incompatible with the version being run (ver',...
               'sion ', tbVer, '). Please re-enter the job specification ',...
               'in the batch window.']);
    end
    
end

%--------------------------------------------------------------------------
% The below functions are currently unused but are useful for comparing
% version numbers. They have been left here only in case they are of use
% for future development. Tom Maullin (09/11/2018)
%--------------------------------------------------------------------------

% Latest compatible version from earliest compatible version hashmap.
function lcv = latComVer(ecv)

    lcv = containers.Map();
    vers = ecv.keys;
    for i = 1:length(vers)
        ver = vers{i};
        
        % If it's an old version we have no recorded version history for
        % earliest or latest compatible versions.
        if isequal(ecv(ver), 'NaN.NaN.NaN')
            lcv(ver) = 'NaN.NaN.NaN';
        % If it's a newer version we need to find the newest version
        % available.
        else
            if ~isKey(lcv, ecv(ver))
                lcv(ecv(ver)) = ver;
            else
                if swe_compareVersions(lcv(ecv(ver)), ver, '<=')
                    lcv(ecv(ver)) = ver;
                end
            end
        end
    end
    
    % Second pass for keys that aren't values in ecv.
    for i = 1:length(vers)
        ver = vers{i};
        
        if ~isKey(lcv, ver)
            lcv(ver) = lcv(ecv(ver));
        end
    end
    
end