function vol=swe_create_vol(fname, DIM, M, varargin)
% FORMAT vol = swe_create_vol(fname, DIM, M [, descrip])
% Initialise a new volume for writing
% 
% fname   - Filename of new image
% DIM     - Row vector giving image dimensions
% M       - 4x4 homogeneous transformation, from V.mat
% descrip - Description to enter into image header
%
%_______________________________________________________________________
% SwE-toolbox

    
if nargin > 3
    descrip = varargin{1};
else
    descrip = '';
end

vol = struct(...
  'fname',    fname,...
  'dim',      DIM',...
  'dt',       [spm_type('float32') spm_platform('bigend')],...
  'mat',      M,...
  'pinfo',    [1 0 0]',...
  'descrip',  descrip);
vol = spm_create_vol(vol);
