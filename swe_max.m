function [N, N_boxcox, Z, M, A, XYZ] = swe_max(X, locationsInVoxels, boxcoxInfo)
  % Sizes, local maxima and locations of excursion sets on volume data
  % FORMAT [N, N_boxcox, Z, M, A, XYZ] = spm_max(X, locationsInVoxels, boxcoxInfo)
  %
  % X                     		- an [nx1] array of stat values
  % locationsInVoxels			    - an [3xn] array of locations {in voxels}
  %	boxcoxInfo								- a struct containing info about the Box-Cox normalisation
  %
  % N                     		- a [px1] size of connected components {in voxels}
  % N_boxcox               	  - a [px1] size of normalised cluster size based on a boxcox transformation normalised by the median (Q2) and the upperHalf-IQR (Q3-Q2)
  % Z                     		- stat values of maxima
  % M                     		- location of maxima {in voxels}
  % A                     		- region number
  % XYZ                   		- cell array of voxel locations {in voxels}
  %
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$

  [N, Z, M, A, XYZ]  = spm_max(X, locationsInVoxels);

  try
    scalingFactorNorm = swe_invNcdf(0.75);
    N_boxcox = (swe_boxCoxTransform(N, boxcoxInfo.volume.lambda) - boxcoxInfo.volume.median) * (scalingFactorNorm / boxcoxInfo.volume.upperHalfIqr);
  catch
    N_boxcox = [];
  end

end
