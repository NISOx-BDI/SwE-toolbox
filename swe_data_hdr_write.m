function V = swe_data_hdr_write(fname, DIM, M, descrip, metadata, varargin)
  % Initialise a new file for writing
  % =========================================================================
  % FORMAT V = swe_data_hdr_write(fname, DIM, M, descrip, metadata[, dataType])
  % -------------------------------------------------------------------------
  % Inputs:
  %   - fname:    Filename of new image
  %   - DIM:      Row vector giving image dimensions
  %   - M:        4x4 homogeneous transformation, from V.mat
  %   - descrip:  Description to enter into image header
  %   - metadata: metadata from GIfTI file (SPM set metadata = {} for NIfTI)
  %   - dataType: data format (e.g., 'float32')
  % =========================================================================
  % Version Info:  $Format:%ci$ $Format:%h$
  if nargin > 5
    dataType = varargin{1};
  else
    dataType = 'float32';
  end

  V = struct(...
    'fname',    fname,...
    'dim',      DIM,...
    'dt',       [spm_type(dataType) spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  descrip,...
    metadata{:});

  if isfield(V, 'ciftiTemplate')
    [~, sliceInd] = swe_get_file_extension(V.ciftiTemplate);
    if isempty(sliceInd)
      sourceName = V.ciftiTemplate;
    else
      sourceName = V.ciftiTemplate(1:( end - numel(sliceInd) ));
    end

    % make sure we select only one slice
    V = swe_data_hdr_read(sprintf('%s,1', sourceName), true);

    V.fname = fname;
    V.private.dat.fname = fname;
    V.private.dat = file_array(fname,...
                                 [1,1,1,1,1,V.dim(1)],...
                                 V.private.dat.dtype,...
                                 0,...
                                 1,...
                                 0);
    V.n = [1 1];
    % modify the xml if the number of point is > 1
    xml = char(V.private.hdr.ext.edata(:)');
    ind = strfind(xml, 'NumberOfSeriesPoints="');
    if ~isempty(ind)
      ind  = ind  + numel('NumberOfSeriesPoints="');
      xml(ind) = '1';
      while ~strcmp(xml(ind+1), '"')
        xml(ind+1) = [];
      end
    end

    % modify the xml if there is more than one <NamedMap>
    ind = strfind(xml, '<NamedMap>');
    if numel(ind > 1)
      ind2 = strfind(xml,'</NamedMap>');
      sliceInd = str2num(sliceInd);
      % remove the <NamedMap>'s after the selected slice unless this is the last
      if sliceInd ~= numel(ind)
        xml(ind(sliceInd+1):(ind2(end) + numel('</NamedMap>'))) = [];
      end
      % remove the <NamedMap>'s before the selected slice unless this is the first
      if sliceInd ~= 1
        xml(ind(1):(ind2(sliceInd-1) + numel('</NamedMap>'))) = [];
      end
    end
    V.private.hdr.ext.edata = uint8(xml)';

    create(V.private)
    V.private.dat(:) = 0;
    V = swe_data_hdr_read(fname, false);
    V.descrip = descrip;
  else
    V = spm_data_hdr_write(V);
  end

end
