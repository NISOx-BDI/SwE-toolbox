% =========================================================================
% Threshold an image that need to be selected by the user from a selection
% window
% =========================================================================
% FORMAT: spm_read_vols(threshold, minimumClusterSize)
% -------------------------------------------------------------------------
% Inputs:
%
%  - threshold - minimum value (inclusive) of the surviving voxels
%  - minimumClusterSize - the minimum size ofthe surviving clusters
% =========================================================================
function swe_thresholdImage(threshold, minimumClusterSize)

  inputImageName = spm_select(1, 'image');
  [pth, bnm, ext] = spm_fileparts(inputImageName);
  VI = spm_vol(inputImageName);
  [Z, XYZ] = spm_read_vols(VI);
  XYZ = inv(VI.mat) * [XYZ; ones(1,VI.dim(1)*VI.dim(2)*VI.dim(3))];
  XYZ = round(XYZ(1:3,:));

  Z = spm_get_data(inputImageName, XYZ);
  VI.fname = fullfile(pth, [bnm '_thresholded' ext]);
  VI.descrip = [VI.descrip sprintf(' thresholdValue: %fminimumClusterSize: %i', threshold, minimumClusterSize)];
  VI = spm_create_vol(VI);

  %-Calculate height threshold filtering
  %--------------------------------------------------------------------------
  Q      = find(Z >= threshold);

  %-Apply height threshold
  %--------------------------------------------------------------------------
  Z      = Z(:,Q);
  XYZ    = XYZ(:,Q);
  if isempty(Q)
    fprintf('\n');                                                      %-#
    warning('SwE:NoVoxels','No voxels survive the thresholding');
  else

    %-Calculate extent threshold filtering
    %----------------------------------------------------------------------
    A     = spm_clusters(XYZ);
    Q     = [];
    for i = 1:max(A)
      j = find(A == i);
      if length(j) >= minimumClusterSize; Q = [Q j]; end
    end

    % ...eliminate voxels
    %----------------------------------------------------------------------
    Z     = Z(:,Q);
    XYZ   = XYZ(:,Q);
    if isempty(Q)
      fprintf('\n');                                                  %-#
      warning('SwE:NoVoxels','No voxels survive the thresholding');
    end

  end

  tmp= nan(VI.dim);
  Q = cumprod([1,VI.dim(1:2)])*XYZ - ...
    sum(cumprod(VI.dim(1:2)));
  tmp(Q) = Z;
  spm_write_vol(VI, tmp);


