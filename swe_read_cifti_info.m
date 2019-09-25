function [surfaces, volume, volumes] = swe_read_cifti_info(filename)
  % Get the file extension
  % =========================================================================
  % FORMAT swe_get_file_extension(filename)
  % -------------------------------------------------------------------------
  % Inputs:
  %   - filename: the name of the cifti file
  % Outputs:
  %   - surfaces: information from each surface brain structure
  %   - volume: information from the merged volumetric brain structure
  %   - volumes: information from each volumetric brain structure
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$
	
	% fetch the cifti information
	ciftiObject  = swe_cifti(filename);
  xml = char(ciftiObject.hdr.ext.edata(:)');
  tree = xmltree(xml);
  root_uid = root(tree);
  hdr = convert(tree);
  hdr.attributes = getAttributes(tree, root_uid);
  
  % Check if input is a CIFTI-2 file
  switch hdr.attributes.Version
    case {'2'}
      % ok
    case {'1','1.0'}
      warning('CIFTI-1 file not supported.');
      return;
    otherwise
      warning('Unknown CIFTI version: %s.', hdr.attributes.Version);
      return;
  end
  
  % extract surfaces and a single volume (by concatenating all brain volumes)
  surfaces = {};
  volumes = {};
  volume = struct;
  volume.XYZ = [];
  volume.indices = [];
  uid_maps = find(tree, 'CIFTI/Matrix/MatrixIndicesMap');
  for i = 1:numel(hdr.Matrix.MatrixIndicesMap)
    map = hdr.Matrix.MatrixIndicesMap{i};
    map.attributes = getAttributes(tree, uid_maps(i));
    switch map.attributes.AppliesToMatrixDimension
      case '0'
        % for now, do nothing
      case '1'
        if ~strcmp(map.attributes.IndicesMapToDataType, 'CIFTI_INDEX_TYPE_BRAIN_MODELS')
          error('Mismatch between NIfTI intent code and CIFTI index type.');
        end
        if isstruct(map.BrainModel), map.BrainModel = {map.BrainModel}; end
        uid_brainModels = find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel');
        for j = 1:numel(map.BrainModel)
          map.BrainModel{j}.attributes = getAttributes(tree, uid_brainModels(j));
          switch map.BrainModel{j}.attributes.ModelType
            case 'CIFTI_MODEL_TYPE_SURFACE'
              nV  = str2double(map.BrainModel{j}.attributes.SurfaceNumberOfVertices);
              off = str2double(map.BrainModel{j}.attributes.IndexOffset);
              iV  = str2num(map.BrainModel{j}.VertexIndices) + 1;
              if numel(iV) ~= str2double(map.BrainModel{j}.attributes.IndexCount)
                error('Problem with number of vertices.');
              end
              BrainStructure = map.BrainModel{j}.attributes.BrainStructure;
              surfaces{end+1} = struct('nV', nV, 'off', off, 'iV', iV, 'brainStructure', BrainStructure);
            case 'CIFTI_MODEL_TYPE_VOXELS'
              off = str2double(map.BrainModel{j}.attributes.IndexOffset);
              iV  = str2num(map.BrainModel{j}.VoxelIndicesIJK) + 1;
              iV  = reshape(iV,3,[]);
              if size(iV, 2) ~= str2double(map.BrainModel{j}.attributes.IndexCount)
                error('Problem with number of voxels.');
              end
							indices = (off + 1):(off + str2num(map.BrainModel{j}.attributes.IndexCount));
							BrainStructure = map.BrainModel{j}.attributes.BrainStructure;
							volumes{end+1} = struct('off', off, 'iV', iV, 'brainStructure', BrainStructure);
              volume.indices = [volume.indices indices];
              volume.XYZ = [volume.XYZ iV];
          end
        end
        if ~isempty(volume)
          uid_volume = find(tree, 'CIFTI/Matrix/MatrixIndicesMap/Volume');
          map.Volume.attributes = getAttributes(tree, uid_volume);
          volume.DIM = str2num(map.Volume.attributes.VolumeDimensions);
          volume.M = str2num(map.Volume.TransformationMatrixVoxelIndicesIJKtoXYZ);
          volume.M = reshape(volume.M, [4 4])' * [eye(4,3) [-1 -1 -1 1]'];
          % assumes MeterExponent is -3
        end
        otherwise
          error('Data have to be dense scalars or series.');
    end
  end
end 
%==========================================================================
function attrb = getAttributes(tree, uid)
    attrb = attributes(tree, 'get', uid);
    if ~isstruct(attrb), attrb = [attrb{:}]; end
    attrb = cell2struct({attrb.val},{attrb.key},2);
end