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

  if isa(V(1).private, 'cifti')
    if ~isstruct(V)
      V = swe_data_hdr_read(V);
    end
    Y = zeros(numel(V), prod(V(1).dim));
    if isempty(varargin)
      for i=1:numel(V)
        Y(i,:) = V(i).private.dat(:);
      end
    else
      indices = varargin;
      if strcmpi(indices{1},'xyz')
        indices = {indices{2}(1,:)};
      end
      n = get_ndata(V(1).dim,indices{:});
      Y = zeros(numel(V),prod(n));
      for i=1:numel(V)
        if numel(indices) == 1
            ind = {indices{1} + (V(i).n(1)-1)*prod(V(i).dim)};
        else
            ind = indices;
        end
        Y(i,:) = reshape(V(i).private.dat(ind{:}),1,[]);
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

%==========================================================================
function n = get_ndata(dim,varargin)
  n = zeros(1,numel(varargin));
  for i=1:numel(varargin)
      if isequal(varargin{i},':')
          if i==numel(varargin)
              n(i) = dim(i); %prod(dim(i:end));
          else
              n(i) = dim(i);
          end
      else
          n(i) = numel(varargin{i});
      end
  end
end