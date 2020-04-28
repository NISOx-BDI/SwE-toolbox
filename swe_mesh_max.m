function [N, N_area, N_boxcox, Z, M, A, XYZ] = swe_mesh_max(X, locationsInVertices, G, boxcoxInfo, areaFile)
  % Sizes, local maxima and locations of excursion sets on volume data
  % FORMAT [N, N_area, N_boxcox, Z, M, A, XYZ] = spm_max(X, locationsInVoxels, G, boxcoxInfo, areaFile)
  %
  % X                     		- an [nx1] array of stat values
  % locationsInVertices			  - an [1xn] array of locations {in vertices}
  %	boxcoxInfo								- a struct containing info about the Box-Cox normalisation
  % areaFile									- optional path to a file containing the area of ech vertex
  %
  % N                     		- a [px1] size of connected components {in vertices}
  % N_area                    - a [px1] surface area size of connected components {in vertices}
  % N_boxcox               	  - a [px1] size of normalised cluster size based on a boxcox transformation normalised by the median (Q2) and the upperHalf-IQR (Q3-Q2)
  % Z                     		- stat values of maxima
  % M                     		- location of maxima {in vertices}
  % A                     		- region number
  % XYZ                   		- cell array of voxel locations {in vertices}
  %
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$

  [N, Z, M, A, XYZ] = spm_mesh_max(X, locationsInVertices, G);
  
  canComputeArea = nargin > 4 && ~isempty(areaFile);
  if canComputeArea
    N_area = zeros(numel(N),1);
    for ii = 1:max(A)
      N_area(A == ii) = sum(swe_data_read(areaFile, XYZ{ii}(1,:)));
    end
  else
    N_area = [];
  end

  canComputeBoxCoxNorm = nargin > 3 && ~isempty(boxcoxInfo);
  if canComputeBoxCoxNorm
    scalingFactorNorm = swe_invNcdf(0.75);
    if canComputeArea
      N_boxcox = (swe_boxCoxTransform(N_area, boxcoxInfo.surfaces.lambda) - boxcoxInfo.surfaces.median) * (scalingFactorNorm / boxcoxInfo.surfaces.upperHalfIqr);
    else
      N_boxcox = (swe_boxCoxTransform(N, boxcoxInfo.surfaces.lambda) - boxcoxInfo.surfaces.median) * (scalingFactorNorm / boxcoxInfo.surfaces.upperHalfIqr);
    end
  else
    N_boxcox = [];
  end

end