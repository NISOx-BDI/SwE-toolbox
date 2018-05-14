function vol=swe_create_vol(fname, dim, m, varargin)
% FORMAT vol = swe_create_vol(fname, dim, m [, desc])
% Initialise a new volume for writing
% 
% fname    - Filename of new image
% dim      - Row vector giving image dimensions
% m        - 4x4 homogeneous transformation, from V.mat
% desc     - Description to enter into image header
% meshData - Boolean stating whether we output gifti or not.
%
%_______________________________________________________________________
% SwE-toolbox
    
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

vol = deal(struct(...
  'fname',    fname,...
  'dim',      dim',...
  'dt',       [spm_type('float32') spm_platform('bigend')],...
  'mat',      m,...
  'pinfo',    [1 0 0]',...
  'descrip',  descrip));

if meshData
    vol = spm_data_hdr_write(vol);
else
    vol = spm_create_vol(vol);
end



