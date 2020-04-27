classdef swe_DataType
  % Class mimicking an enumeration class with values:
	% Nifti, Gifti, Cifti, Mat, VolumeMat, SurfaceMat.
  %
  % It is used instead of a Matlab enumeration class due to the lack of 
  % enumeration support in octave
  %
  % FORMAT dataType = swe_DataType(dataTypeChar)
  % dataTypeChar - One of the following char arrays:
  %              'Nifti', 'Gifti', 'Cifti', 'Mat', 'VolumeMat','SurfaceMat'
  % dataType     - an object of class swe_DataType
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$
  
  properties (Constant)
    enum = {'Nifti', 'Gifti', 'Cifti', 'Mat', 'VolumeMat', 'SurfaceMat'};
  end
  
  properties 
    value
  end
  
  methods
    
    function obj = swe_DataType(dataType)
      if nargin == 0
        error( "swe_DataType cannot be initialized empty; choose from 'Nifti', 'Gifti', 'Cifti', 'Mat', 'VolumeMat' or 'SurfaceMat'.");
      end
      if ~ismember(dataType, obj.enum)
        error("%s is not a valid swe_DataType. Choose from 'Nifti', 'Gifti', 'Cifti', 'Mat', 'VolumeMat' or 'SurfaceMat'.", dataType);
      end
      obj.value = dataType;
    end
    
    function display(obj)
      display(get(obj));
    end

    function out = get(obj)
      out = obj.value;
    end

    function out = eq(a,b)
      out = strcmp(a.value, b.value);
    end

  end

end