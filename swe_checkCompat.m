% matVer - version of SwE recorded in the `.mat` file.
% tbVer  - version of SwE toolbox currently being run.
function swe_checkCompat(matVer, tbVer)

    if isequal(swe('ver'),SwE.ver)
       return; 
    end
    
    % The below hashmap records the earliest compatable version for each
    % release of the SwE toolbox. I.e. when making a new release, say you
    % are releasing version "y.y.y", please set earliestCompatVer("y.y.y")
    % equal to "x.x.x" where "x.x.x" is the oldest version of the toolbox
    % which "y.y.y" can accept `SwE.mat` files from.
    earliestCompatVer = containers.Map();
    
    % These versions did not record version numbers in the SwE.mat file and
    % therefore cannot be checked in the same way. However, none of these
    % should be compatabile with version 2.0.1 (the version in which this
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
    
    % Record earliest compatable versions.
    earliestCompatVer('2.0.0') = '2.0.0';
    earliestCompatVer('2.0.1') = '2.0.1';

end