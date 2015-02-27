function vechCovVechV = swe_vechCovVechV(covVis,dofMat,type)
% Compute vech(Cov(vech(\hat V))) following one of the 3 following
% methodologies:
% type = 1: based on the estimate "II" proposed in Guillaume (2015).
%           It accounts partially for missing data (not for a missing data 
%           bias) and for a small sample bias
% type = 2: based on the estimate "III" proposed in Guillaume (2015).
%           It accounts for missing data (included the missing data bias), 
%           but not for the small sample bias.
% type = 3: based on the estimate "II" proposed in Guillaume (2015), but
%           requires no missing data (so accounts for the small sample bias).
% covVis: matrix (nCovVis X nVox)  containing vech(\hat V) for several voxels 
% dofMat: for type 1 and 2, matrix (nCovVis X nCovVis) containing degrees of freedom 
%         information (precalculated in swe_cp.m). 
%         for type 3, scalar.
%         See Guillaume (2015)for more information.
% Note: there are probably ways to compute this quicker, but will do for
% now.
% By Bryan Guillaume
[nCovVis, nVox] = size(covVis);
nVis = -0.5 + sqrt(0.25 + 2 * nCovVis);
nVechCovVechV      = nCovVis * (nCovVis + 1) / 2;
vechCovVechV       = zeros(nVechCovVechV,nVox);
[cIndex, rIndex] = find(tril(ones(nVis))==1);

if type == 3
    factor = dofMat / (dofMat - 1) / (2* dofMat + 1);
end

it = 0;        
for i = 1:nCovVis
    for ii = i:nCovVis
        it = it + 1;
        
        ind11 = rIndex == rIndex(i) & cIndex == rIndex(ii);
        ind12 = rIndex == rIndex(i) & cIndex == cIndex(ii);
        ind22 = (rIndex == cIndex(i) & cIndex == cIndex(ii)) | (cIndex == cIndex(i) & rIndex == cIndex(ii));
        ind21 = (rIndex == cIndex(i) & cIndex == rIndex(ii)) | (cIndex == cIndex(i) & rIndex == rIndex(ii));
        indi1 = rIndex == rIndex(i) & cIndex == rIndex(i);
        indi2 = rIndex == cIndex(i) & cIndex == cIndex(i);
        indii1 = rIndex == rIndex(ii) & cIndex == rIndex(ii); 
        indii2 = rIndex == cIndex(ii) & cIndex == cIndex(ii);
        
        switch type
            
            case 1
                tmp = 1 + 2 * dofMat(i,ii) * dofMat(ind11,ind22) * dofMat(ind12,ind21) - ...
                    dofMat(i,ii) * dofMat(ind11,ind22) - ...
                    dofMat(ind11,ind22) * dofMat(ind12,ind21) - ...
                    dofMat(i,ii) * dofMat(ind12,ind21);
     
                vechCovVechV(it,:) = dofMat(i,ii) / tmp * (...
                    (2 * dofMat(ind11,ind22) * dofMat(ind12,ind21) - dofMat(ind11,ind22) - dofMat(ind12,ind21)) * covVis(i,:) .* covVis(ii,:) + ...
                    (1 - dofMat(ind12,ind21)) * covVis(ind11,:) .* covVis(ind22,:) + ...
                    (1 - dofMat(ind11,ind22)) * covVis(ind12,:) .* covVis(ind21,:));

            case 2
                vechCovVechV(it,:) = dofMat(i,ii) * (covVis(ind11,:).* covVis(ind22,:) + covVis(ind12,:).* covVis(ind21,:)) + ...
                    covVis(i,:) .*( covVis(ind11,:) .* covVis(ind12,:) ./ covVis(indi1,:) * (dofMat(indi1,ii) - dofMat(i,ii)) + ...
                    covVis(ind21,:) .* covVis(ind22,:) ./ covVis(indi2,:) * (dofMat(indi2,ii) - dofMat(i,ii))) + ...
                    covVis(ii,:) .*( covVis(ind11,:) .* covVis(ind21,:) ./ covVis(indii1,:) .* (dofMat(i,indii1) - dofMat(i,ii)) + ...
                    covVis(ind12,:) .* covVis(ind22,:) ./ covVis(indii2,:) * (dofMat(i,indii2) - dofMat(i,ii))) + ...
                    0.5 * covVis(i,:).*covVis(ii,:).*(...
                    covVis(ind11,:).^2./covVis(indi1,:)./covVis(indii1,:) * (dofMat(indi1,indii1) + dofMat(i,ii) - dofMat(i,indii1) - dofMat(indi1,ii)) + ...
                    covVis(ind21,:).^2./covVis(indi2,:)./covVis(indii1,:) * (dofMat(indi2,indii1) + dofMat(i,ii) - dofMat(i,indii1) - dofMat(indi2,ii)) + ...
                    covVis(ind12,:).^2./covVis(indi1,:)./covVis(indii2,:) * (dofMat(indi1,indii2) + dofMat(i,ii) - dofMat(i,indii2) - dofMat(indi1,ii)) + ...
                    covVis(ind22,:).^2./covVis(indi2,:)./covVis(indii2,:) * (dofMat(indi2,indii2) + dofMat(i,ii) - dofMat(i,indii2) - dofMat(indi2,ii)));
            case 3
                vechCovVechV(it,:) = factor *...
                    (2 * dofMat * covVis(i,:) .* covVis(ii,:) -...
                    covVis(ind11,:) .* covVis(ind22,:) - ...
                    covVis(ind12,:) .* covVis(ind21,:));
        end
    end
end
