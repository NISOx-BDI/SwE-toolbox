function obj = swe_DataType(dataType)
  % Create a swe_DataType object
  % Its goal is to mimick the behaviour of an enumaration of data types that is
  % compatible with Octave.
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$

  availableDataTypes = {'Nifti', 'Gifti', 'Cifti', 'Mat', 'VolumeMat', 'SurfaceMat'};

  if nargin == 0
    error( "swe_DataType cannot be initialized empty; choose from 'Nifti', 'Gifti', 'Cifti', 'Mat', 'VolumeMat' or 'SurfaceMat'.");
  end

  if ~ismember(dataType, availableDataTypes)
    error("%s is not a valid swe_DataType. Choose from 'Nifti', 'Gifti', 'Cifti', 'Mat', 'VolumeMat' or 'SurfaceMat'.", dataType);
  end

  obj.value = dataType;

  obj = class(obj, 'swe_DataType');

end
