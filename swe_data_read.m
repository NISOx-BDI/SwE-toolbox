function Y = swe_data_read(V,varargin)
  % Read data from disk [Y = V(I)]
  % FORMAT Y = swe_data_read(V)
  % V        - a structure array (see swe_data_hdr_read)
  % Y        - an array of data values; the last dimension indexes numel(V)
  %
  % FORMAT Y = swe_data_read(V,'slice',S)
  % V        - a structure array of image volumes (see swe_data_hdr_read)
  % S        - an array of slice indices
  % Y        - an array of data values with dimensions (x,y,s,v)
  %
  % FORMAT Y = swe_data_read(V,'xyz',XYZ)
  % V        - a structure array (see swe_data_hdr_read)
  % XYZ      - a [n x m] array of m coordinates {voxel (n=3 or 4)/vertex (n=1)}
  % Y        - an array of data values with dimensions (v,m)
  %
  % FORMAT Y = swe_data_read(V,I1,I2,...)
  % V        - a structure array (see swe_data_hdr_read)
  % I1,I2,...- subscript arrays
  % Y        - an array of data values with dimensions (v,m)
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$

  if ~isstruct(V)
  	V = swe_data_hdr_read(V);
  end
    
  if isa(V(1).private, 'swe_cifti')
    if isempty(varargin)
      Y = zeros(prod(V(1).dim), numel(V)); % to be coherent with spm_read_vols and gifti format
      for i=1:numel(V)
        Y(:,i) = V(i).private.dat(1,1,1,1,V(i).n(1),:);
      end
    else
      indices = varargin;
      if strcmpi(indices{1},'xyz')
        indices = {indices{2}(1,:)};
      end
      Y = zeros(numel(V),numel(indices{:}));
      for i=1:numel(V)
        Y(i,:) = V(i).private.dat(1,1,1,1,V(i).n(1),indices{:});
      end
    end
  else
    if isempty(varargin)
      Y = spm_data_read(V);
    else
      Y = spm_data_read(V,varargin{:});
    end
  end
end
