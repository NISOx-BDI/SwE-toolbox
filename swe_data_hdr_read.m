function V = swe_data_hdr_read(P, varargin)
  % Read data information from file(s)
  % FORMAT V = spm_data_hdr_read(P)
  % P        - a char or cell array of filenames
  % V        - a structure array containing data information
  %
  % This function behaves like spm_data_hdr_read but can also read the
  % headers of CIfTI files.
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$
  P2 = cellstr(P);
  file_ext = swe_get_file_extension(P2{1});
  isCifti  = strcmpi(file_ext,'.dtseries.nii') ||  strcmpi(file_ext,'.dscalar.nii');

  if nargin > 1
    readExt = varargin{1};
  else
    readExt = false;
  end

  if isCifti
    it = 1;
    for i=1:numel(P2)
      [file_ext, sliceInd] = swe_get_file_extension(P2{i});
      if ~(strcmpi(file_ext,'.dtseries.nii') ||  strcmpi(file_ext,'.dscalar.nii'))
        error('Not all files are CIfTI files');
      end

      % if some slices of data are requested, remove this from the filename
      if ~isempty(sliceInd)
        P2{i} = P2{i}(1:(end-numel(sliceInd)));
      end

      ciftiObject = swe_cifti(P2{i}, readExt);

      if isempty(sliceInd)
        sliceInd = 1:ciftiObject.dat.dim(5);
      else
        sliceInd = str2num(sliceInd);
        if any(sliceInd > ciftiObject.dat.dim(5))
          error('At least one slice index exceeds the dimension of the data array!');
        end
      end

      for iSlice = sliceInd
        V(it) = default_hdr_struct;
        V(it).fname    = P2{i};
        V(it).private  = ciftiObject;
        V(it).dim      = [ciftiObject.dat.dim(6), 1];
        V(it).n = [iSlice 1];
        it = it + 1;
      end
    end
  else
    V = spm_data_hdr_read(P);
  end
end

function V = default_hdr_struct
	V = struct(...
			'fname',   '',...
			'dim',     [0 0 0],...
			'dt',      [spm_type('float64') spm_platform('bigend')],...
			'pinfo',   [1 0 0]',...
			'mat',     eye(4),...
			'n',       [1 1],...
			'descrip', '',...
			'private', []);
end
