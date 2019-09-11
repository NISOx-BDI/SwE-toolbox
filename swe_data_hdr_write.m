function V = swe_data_hdr_write(fname, DIM, M, descrip, metadata, varargin)
  % Initialise a new file for writing
  % =========================================================================
  % FORMAT V = swe_data_hdr_write(fname, DIM, M, descrip, metadata[, dataType])
  % -------------------------------------------------------------------------
  % Inputs: 
  %   - fname:    Filename of new image
  %   - DIM:      Row vector giving image dimensions
  %   - M:        4x4 homogeneous transformation, from V.mat
  %   - descrip:  Description to enter into image header
  %   - metadata: metadata from GIfTI file (SPM set metadata = {} for NIfTI)
  %   - dataType: data format (e.g., 'float32')
  % =========================================================================
  % Version Info:  $Format:%ci$ $Format:%h$
  if nargin > 5
    dataType = varargin{1};
  else
    dataType = 'float32';
  end

  V = struct(...
    'fname',    fname,...
    'dim',      DIM,...
    'dt',       [spm_type(dataType) spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  descrip,...
    metadata{:});

  V = spm_data_hdr_write(V);
      
end
  