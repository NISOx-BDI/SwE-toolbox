function fileExtension = swe_get_file_extension(filename, varargin)
  % Get the file extension
  % =========================================================================
  % FORMAT swe_get_file_extension(filename)
  % -------------------------------------------------------------------------
  % Inputs:
  %   - filename: the name of the file
  % Outputs:
  %   - fileExtension: extension of the file taken after the first '.' found
  % =========================================================================
  % FORMAT swe_get_file_extension(filename, lastDot)
  % -------------------------------------------------------------------------
  % Inputs:
  %   - filename: the name of the file
  %   - lastDot: if true, the extension is taken after the last '.'
  %              if false, the extension is taken after the first '.'
  % Outputs:
  %   - fileExtension: extension of the file
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$
  
  if nargin == 1
    lastDot = false;
  elseif nargin == 2
    lastDot = varargin{1};
  else
    error('too many inputs');
  end
  
  indexDots = find(filename == '.');

  if lastDot
    fileExtension = filename(indexDots(end):end);
  else
    fileExtension = filename(indexDots(1):end);
  end

end