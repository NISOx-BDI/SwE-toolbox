function V = swe_data_hdr_read(P)
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
  P = cellstr(P);
  file_ext = swe_get_file_extension(P{1});
  isCifti  = strcmpi(file_ext,'.dtseries.nii') ||  strcmpi(file_ext,'.dtscalar.nii');
  if isCifti
    [V(1:numel(P),1)] = deal(default_hdr_struct);
    for i=1:numel(P)
      V(i).fname    = P{i};
      V(i).private  = swe_cifti(P{i}, false);
      V(i).dim      = size(squeeze(V(i).private.dat));
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