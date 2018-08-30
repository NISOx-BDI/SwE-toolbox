% =========================================================================
% Initialise a new volume for writing
% =========================================================================
% FORMAT vol = swe_create_vol(fname, DIM, M [, descrip])
% -------------------------------------------------------------------------
% Inputs: 
%   - fname:    Filename of new image
%   - DIM:      Row vector giving image dimensions
%   - M:        4x4 homogeneous transformation, from V.mat
%   - descrip:  Description to enter into image header
%   - meshData: Boolean stating whether we output gifti or not.
% =========================================================================
% SwE-toolbox

function vol=swe_create_vol(fname, DIM, M, varargin)
    
if nargin > 3
    descrip = varargin{1};
else
    descrip = '';
end

if nargin > 4
    meshData = varargin{2};
else
    meshData = false;
end

vol = struct(...
  'fname',    fname,...
  'dim',      DIM',...
  'dt',       [spm_type('float32') spm_platform('bigend')],...
  'mat',      M,...
  'pinfo',    [1 0 0]',...
  'descrip',  descrip);

if meshData
    vol = spm_data_hdr_write(vol);
else
    vol = spm_create_vol(vol);
end
