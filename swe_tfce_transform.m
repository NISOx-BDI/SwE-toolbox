function [tfced] = swe_tfce_transform(img,H,E,C,dh)
% This function calculates the TFCE statistical from a Z statistic image.
% =========================================================================
%   FORMAT [tfced] = matlab_tfce_transform(img,H,E,C,ndh) 
% -------------------------------------------------------------------------
% Inputs:
%   - img the 3D image to be transformed
%   - H height exponent
%   - E extent exponent
%   - C connectivity
%   - dh size of steps for cluster formation
% =========================================================================
% Performs threshold free cluster enhancement on 'img' as per Smith &
% Nichols (2009).
%
% This function was modified from 'matlab_tfce_transform' from MatlabTFCE.

% set cluster thresholds
threshs = 0:dh:max(img(:));
threshs = threshs(2:end);
ndh = length(threshs);

% find positive voxels (greater than first threshold)
nvox = length(img(:));

% find connected components
vals = zeros(nvox,1);
cc = arrayfun(@(x) bwconncomp(bsxfun(@ge,img,x),C), threshs);
for h = 1:ndh
    clustsize = zeros(nvox,1);
    ccc = cc(h);
    voxpercc = cellfun(@numel,ccc.PixelIdxList);
    for c = 1:ccc.NumObjects
        clustsize(ccc.PixelIdxList{c}) = voxpercc(c);
    end
    % calculate transform
    curvals = (clustsize.^E).*(threshs(h)^H);
    vals = vals + curvals;
end
tfced = NaN(size(img));
tfced(:) = vals.*dh;

end

