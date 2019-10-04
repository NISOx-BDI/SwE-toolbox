function [A, clusterSizesInSurfaces, clusterSizesInVolume] = swe_cifti_clusters(ciftiInformation, indSurvivingInCifti)
  % Compute the clusters for CIfTI data
  % FORMAT A = swe_cifti_clusters(ciftiInformation, indSurvivingInCifti)
  % ciftiInformation      	- cifti information
  % indSurvivingInCifti     - an array of locations {in cifti indices}
  % A                       - a vector of cluster indices
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$
  A = zeros(1, numel(indSurvivingInCifti));
  clusterSizesInSurfaces = [];
  if numel(ciftiInformation.surfaces) > 0
    for i = 1:numel(ciftiInformation.surfaces)
      % work out the position of the surviving vertices for this surface
      indInSurface = ciftiInformation.surfaces{i}.off + (1:numel(ciftiInformation.surfaces{i}.iV));
      [isSurviving, indInA] = ismember(indInSurface, indSurvivingInCifti);
      indInA = indInA(isSurviving);
      indSurvivingInSurface = ciftiInformation.surfaces{i}.iV(isSurviving);
      T = false(ciftiInformation.surfaces{i}.nV, 1);
      T(indSurvivingInSurface) = true;
      G = export(gifti(ciftiInformation.surfaces{i}.geomFile), 'patch');
      tmp = spm_mesh_clusters(G,T)';
      tmp = tmp(indSurvivingInSurface);
      A(indInA) = tmp + max(A);
      nClusters = max(tmp);
      if isfield(ciftiInformation.surfaces{i}, 'areaFile')
        area = swe_data_read(ciftiInformation.surfaces{i}.areaFile, indSurvivingInSurface);
        clusterSizesInThisSurface = zeros(1, nClusters);
        for iCluster = 1:nClusters
          clusterSizesInThisSurface(iCluster) = sum(area(tmp == iCluster));
        end
        clusterSizesInSurfaces = [clusterSizesInSurfaces, clusterSizesInThisSurface];
      else
        clusterSizesInSurfaces = [clusterSizesInSurfaces, histc(tmp,1:nClusters)];
      end
    end
  end
  % deal with volumetric data
  if numel(ciftiInformation.volume) > 0
    % work out the position of the surviving voxels
    [isSurviving, indInA] = ismember(ciftiInformation.volume.indices, indSurvivingInCifti);
    indInA = indInA(isSurviving);
    inMask_vol_XYZ = ciftiInformation.volume.XYZ(:, isSurviving);
    tmp = spm_clusters(inMask_vol_XYZ);
    A(indInA) = tmp + max(A);
    nClusters = max(tmp);
    clusterSizesInVolume = histc(tmp,1:nClusters);
  else
    clusterSizesInVolume = [];
  end
end