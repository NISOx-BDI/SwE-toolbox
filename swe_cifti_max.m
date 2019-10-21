function [N, N_area, N_boxcox1, N_boxcox2, Z, M, A, XYZ, Mmm, brainStructureShortLabels, brainStructureLongLabels] = swe_cifti_max(X, indSurvivingInCifti, ciftiInformation, boxcoxInfo)
  % Sizes, local maxima and locations of excursion sets on a cifti file
  % FORMAT [N, N_area, N_boxcox1, N_boxcox2, Z, M, A, XYZ, Mmm] = spm_mesh_max(X, indSurvivingInCifti, ciftiInformation, boxcoxInfo)
  % X                     		- an [nx1] array of stat values
  % indSurvivingInCifti   		- an [nx1] array of locations {in cifti indices}
  % ciftiInformation      		- cifti information
  %
  % N                     		- a [px1] size of connected components {in voxels/vertices}
  % N_area                		- a [px1] size of cluster area {in mm^2} (empty for volume)
  % N_boxcox1               	- a [px1] size of normalised cluster size based on a boxcox transformation normalised by the mean and the std 
  % N_boxcox2               	- a [px1] size of normalised cluster size based on a boxcox transformation normalised by the median and twice (Q3-Q2)
  % Z                     		- stat values of maxima
  % M                     		- location of maxima {in vertices or voxels}
  % A                     		- region number
  % XYZ                   		- cell array of vertices locations in voxels/vertices
  % Mmm                   		- location of maxima {in vertices or voxels}
  % brainStructureShortLabels - short brain structure labels of each maxima
  % brainStructureLongLabels  - long brain structure labels of each maxima
  %__________________________________________________________________________
  N = []; N_area = []; N_boxcox1 = []; N_boxcox2 = [];
  Z = []; M = []; A = []; XYZ = []; Mmm = []; 
  brainStructureShortLabels = []; brainStructureLongLabels = [];
  
  if numel(X) ~= numel(indSurvivingInCifti)
    error('X and indSurvivingInCifti does not contain the same number of elements!')
  end

  if isempty(indSurvivingInCifti)
    return;
  end

  %-Ensure that indSurvivingInCifti contains exactly integers
  %--------------------------------------------------------------------------
  indSurvivingInCifti = round(indSurvivingInCifti);

  %-Detect orientation of X
  %--------------------------------------------------------------------------
  if size(X,1) > 1
    isXColumnVector = true;
  else
    isXColumnVector = false;
  end

  maxA = 0;
  if numel(ciftiInformation.surfaces) > 0
    for i = 1:numel(ciftiInformation.surfaces)
      % work out the position of the surviving vertices for this surface
      indInSurface = ciftiInformation.surfaces{i}.off + (1:numel(ciftiInformation.surfaces{i}.iV));
      [isSurviving, indInX] = ismember(indInSurface, indSurvivingInCifti);
      indInX = indInX(isSurviving);
      indSurvivingInSurface = ciftiInformation.surfaces{i}.iV(isSurviving);
      G = export(gifti(ciftiInformation.surfaces{i}.geomFile), 'patch');
      [N_tmp, Z_tmp, M_tmp, A_tmp, XYZ_tmp] = spm_mesh_max(X(indInX), indSurvivingInSurface, G);
      if ~isempty(A_tmp)
        N = [N; N_tmp];
        if isXColumnVector
          Z = [Z; Z_tmp];
        else
          Z = [Z, Z_tmp'];
        end
        % if there is area info, compute the areas and do the boxcox on the
        % areas. Otherwise, do the boxcox on the number of vertices
        if isfield(ciftiInformation.surfaces{i}, 'areaFile')
          N_area_tmp = zeros(numel(N_tmp),1);
          for ii = 1:max(A_tmp)
            N_area_tmp(A_tmp == ii) = sum(swe_data_read(ciftiInformation.surfaces{i}.areaFile, XYZ_tmp{ii}(1,:)));
          end
          N_area = [N_area; N_area_tmp];
          if ~isempty(boxcoxInfo)
            tmp = boxcox(boxcoxInfo.surfaces.lambda, N_area_tmp);
          end
        else
          if ~isempty(boxcoxInfo)
            tmp = boxcox(boxcoxInfo.surfaces.lambda, N_tmp);
          end
        end
        if ~isempty(boxcoxInfo)
          N_boxcox1 = [N_boxcox1; (tmp - boxcoxInfo.surfaces.mean) ./ boxcoxInfo.surfaces.std];
          N_boxcox2 = [N_boxcox2; (tmp - boxcoxInfo.surfaces.median) ./ boxcoxInfo.surfaces.hiqr];
        end
        % need to convert the surface coordinates into CIfTI coordinates
        isMax = ismember(ciftiInformation.surfaces{i}.iV, M_tmp(1,:));
        M = [M, [indInSurface(isMax); ones(2,sum(isMax))]];
        A = [A; A_tmp + maxA];
        for ii = 1:max(A_tmp)
          isInCluster = ismember(ciftiInformation.surfaces{i}.iV, XYZ_tmp{ii}(1,:));
          XYZ_tmp{ii} = [indInSurface(isInCluster); ones(2,sum(isInCluster))];
        end
        XYZ = [XYZ, XYZ_tmp];
        Mmm = [Mmm, G.vertices(M_tmp(1,:),:)'];
        tmpCell = cell(size(A_tmp));
        [tmpCell{:}] = deal(ciftiInformation.surfaces{i}.brainStructure);
        brainStructureLongLabels = [brainStructureLongLabels; tmpCell];
        [tmpCell{:}] = deal(sprintf('S%i', i));
        brainStructureShortLabels = [brainStructureShortLabels; tmpCell];
        maxA = max(A);
      end
    end
  end
  % deal with volumetric data
  if numel(ciftiInformation.volume) > 0
    [isSurviving, indInX] = ismember(ciftiInformation.volume.indices, indSurvivingInCifti);
    indInX = indInX(isSurviving);
    inMask_vol_XYZ = ciftiInformation.volume.XYZ(:, isSurviving);
    [N_tmp, Z_tmp, M_tmp, A_tmp, XYZ_tmp]  = spm_max(X(indInX), inMask_vol_XYZ);
    if ~isempty(A_tmp)
      N = [N; N_tmp];
      if isXColumnVector
        Z = [Z; Z_tmp];
      else
        Z = [Z, Z_tmp];
      end
			% need to convert the volume coordinates into CIfTI coordinates
      isMax = ismember(ciftiInformation.volume.XYZ', M_tmp', 'rows')'; 
      M = [M, [ciftiInformation.volume.indices(isMax); ones(2,sum(isMax))]];
      A = [A; A_tmp + maxA];
      % need to convert the volume coordinates into CIfTI coordinates
      for i = 1:max(A_tmp)
        isInCluster = ismember(ciftiInformation.volume.XYZ', XYZ_tmp{i}', 'rows')';
        XYZ_tmp{i} = [ciftiInformation.volume.indices(isInCluster); ones(2,sum(isInCluster))];
      end
      XYZ = [XYZ, XYZ_tmp];
      Mmm = [ Mmm, ciftiInformation.volume.M(1:3,:) * [M_tmp; ones(1,size(M_tmp,2))] ];
      tmpCell = cell(size(A_tmp));
      [tmpCell{:}] = deal('VOLUME');
      brainStructureLongLabels = [brainStructureLongLabels; tmpCell];
      [tmpCell{:}] = deal('V');
      brainStructureShortLabels = [brainStructureShortLabels; tmpCell];
      if ~isempty(boxcoxInfo)
        tmp = boxcox(boxcoxInfo.volume.lambda, N_tmp);
        N_boxcox1 = [N_boxcox1; (tmp - boxcoxInfo.volume.mean) ./ boxcoxInfo.volume.std];
        N_boxcox2 = [N_boxcox2; (tmp - boxcoxInfo.volume.median) ./ boxcoxInfo.volume.hiqr];
      end
    end
  end
end