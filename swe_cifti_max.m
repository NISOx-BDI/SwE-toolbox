function [N, Z, M, A, XYZ, Mmm, brainStructureShortLabels, brainStructureLongLabels] = swe_cifti_max(X, indSurvivingInCifti, ciftiInformation)
  % Sizes, local maxima and locations of excursion sets on a cifti file
  % FORMAT [N, Z, M, A, XYZ, Mmm] = spm_mesh_max(X, indSurvivingInCifti, ciftiInformation)
  % X                     		- an [nx1] array of stat values
  % indSurvivingInCifti   		- an [nx1] array of locations {in cifti indices}
  % ciftiInformation      		- cifti information
  %
  % N                     		- a [px1] size of connected components {in vertices}
  % Z                     		- stat values of maxima
  % M                     		- location of maxima {in vertices or voxels}
  % A                     		- region number
  % XYZ                   		- cell array of vertices locations in voxels/vertices
  % Mmm                   		- location of maxima {in vertices or voxels}
  % brainStructureShortLabels - short brain structure labels of each maxima
  % brainStructureLongLabels  - long brain structure labels of each maxima
  %__________________________________________________________________________
  N = []; Z = []; M = []; A = []; XYZ = []; Mmm = []; 
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
        M = [M, M_tmp];
        A = [A; A_tmp + maxA];
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
      M = [M, M_tmp];
      A = [A; A_tmp + maxA];
      XYZ = [XYZ, XYZ_tmp];
      Mmm = [ Mmm, ciftiInformation.volume.M(1:3,:) * [M_tmp; ones(1,size(M_tmp,2))] ];
      tmpCell = cell(size(A_tmp));
      [tmpCell{:}] = deal('VOLUME');
			brainStructureLongLabels = [brainStructureLongLabels; tmpCell];
      [tmpCell{:}] = deal('V');
			brainStructureShortLabels = [brainStructureShortLabels; tmpCell];
    end
  end
end