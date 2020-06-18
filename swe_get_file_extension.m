function [fileExtension, sliceInd] = swe_get_file_extension(filename, varargin)
  % Get the file extension
  % =========================================================================
  % FORMAT swe_get_file_extension(filename)
  % -------------------------------------------------------------------------
  % Inputs:
  %   - filename: the name of the file
  % Outputs:
  %   - fileExtension: extension of the file taken after the first '.' found
  %   - sliceInd: comma separated list of slice indices
  % =========================================================================
  % FORMAT swe_get_file_extension(filename, lastDot)
  % -------------------------------------------------------------------------
  % Inputs:
  %   - filename: the name of the file
  %   - lastDot: if true, the extension is taken after the last '.'
  %              if false, the extension is taken after the first '.'
  % Outputs:
  %   - fileExtension: extension of the file
  %   - sliceInd: comma separated list of slice indices
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
  nDots = numel(indexDots);

  if nDots == 0
    error('No extension found!');
  else
    fileExtension = filename(indexDots(end):end);
    % look for a comma in the extension
    ind = find(fileExtension == ',');
    if ~isempty(ind)
      sliceInd = fileExtension(ind(1):end);
      fileExtension = fileExtension(1:(ind(1)-1));
      filename = filename( 1:( end - numel(sliceInd) ) );
    else
      sliceInd = '';
    end
  end

  if nDots > 1 && ~lastDot && strcmpi(fileExtension, '.nii')
    tmp = filename(indexDots(end-1):end);
    if strcmpi(tmp, '.dtseries.nii') || strcmpi(tmp, '.dscalar.nii')
      fileExtension = tmp;
    end
  end
end
