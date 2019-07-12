function result = swe_compareVersions(ver1, ver2, comparisonOperator)
% Compare two version number of the type 'x.x.x'
% =========================================================================
% FORMAT: result = swe_compareVersions(ver1, ver2, comparisonOperator)
% -------------------------------------------------------------------------
% Inputs: 
%   - ver1: first version number of the type 'x.x.x'
%   - ver2: second version number of the type 'x.x.x'
%   - comparisonOperator: comparison operator ('<', '<=', '==', '>=' or '>')
% =========================================================================
% Author: Bryan Guillaume & Tom Maullin (05/07/2019)
% Version Info:  $Format:%ci$ $Format:%h$
  
  % Split version number into parts.
  ver1 = strsplit(ver1, '.');
  ver2 = strsplit(ver2, '.');
  
  % Work out base for comparison.
  base = max(cellfun(@(a) str2num(a), {ver1{:} ver2{:}})) + 1;
  
  % Convert to num.
  ver1 = cellfun(@(a) str2num(a), {ver1{:}});
  ver2 = cellfun(@(a) str2num(a), {ver2{:}});
  
  % Make into number in base `base`
  ver1 = ver1(1)*base^2 + ver1(2)*base^1 + ver1(3)*base^0;
  ver2 = ver2(1)*base^2 + ver2(2)*base^1 + ver2(3)*base^0;
  
  if strcmp(comparisonOperator, '<=')
    result = (ver1 <= ver2);
  elseif strcmp(comparisonOperator, '<')
    result = (ver1 < ver2);
  elseif strcmp(comparisonOperator, '==')
    result = (ver1 == ver2);
  elseif strcmp(comparisonOperator, '>')
    result = (ver1 > ver2);
  elseif strcmp(comparisonOperator, '>=')
    result = (ver1 >= ver2);
  else
    error('unknown comparison operator')
  end
end