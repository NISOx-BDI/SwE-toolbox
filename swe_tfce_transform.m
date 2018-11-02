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
% -------------------------------------------------------------------------
% Additional Octave compatability using SPM functions was added by 
% Tom Maullin (19/09/2018)
%
% Note: Using the SPM functions is almost 10 times as slow as bwconncomp
% and only recommended when the latter is not available.

% check if we are in octave.
bwfuncexists = exist('bwconncomp','file')==2;

% set cluster thresholds
threshs = 0:dh:max(img(:));
threshs = threshs(2:end);
ndh = length(threshs);

% find positive voxels (greater than first threshold)
nvox = length(img(:));

if bwfuncexists

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

% In Octave bwconncomp is not available unless the 'image' package has been
% installed. The same computation can be performed with SPM functions but
% in much slower time.
else
    
    % find connected components
    vals = zeros(size(img));
    for h = 1:ndh

        % Get cluster labels map.
        clustsize = zeros(size(img));
        [clabels, nc]=spm_bwlabel(bsxfun(@ge,img,threshs(h)).*1, C);

        % Find out how many voxels per cluster
        voxpercc = arrayfun(@(x) numel(clabels(clabels==x)), 1:nc);

        % Re-label clusters with their size.
        for c = 1:nc
            clustsize(clabels==c)=voxpercc(c);
        end

        % calculate transform
        curvals = (clustsize.^E).*(threshs(h)^H);
        vals = vals + curvals;
    end
    tfced = vals.*dh;
end

end

