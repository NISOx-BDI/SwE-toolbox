function V = swe_data_write(V, Y, varargin)
  % Write data to disk [V(I) = Y]
  % FORMAT V = swe_data_write(V,Y)
  % V        - a structure array (see swe_data_hdr_read)
  % Y        - an array of data values
  %
  % FORMAT V = swe_data_write(V,Y,I)
  % V        - a structure array (see swe_data_hdr_read)
  % Y        - an array of data values
  % I        - linear index to data values
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$
  
  if isfield(V,'private')
    cl = class(V.private);
  elseif isfield(V,'dat')
    cl = 'struct';
  else
    error('Unkwown data type.');
  end

  if strcmpi(cl, 'cifti')
    if isempty(varargin)
      V.private.dat = subsasgn(V.private.dat, substruct('()',{':'}), Y);
    else
      try
        V.private.dat(varargin{1}) = reshape(Y,size(varargin{1}));
      catch
        V.private.dat(varargin{1}) = reshape(Y,size(varargin{1}))';
      end
    end
  else
    V = spm_data_write(V, Y, varargin);
  end
end