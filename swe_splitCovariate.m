% =========================================================================
% This function splits a time-varying covariate into a cross-sectional
% covariate and a longitudinal covariate. The split is done first by 
% extracting the within-subect average to build the cross-sectional 
% covariate, which is also centered. Second, the longitudinal covariate is 
% simply the difference between the original covariate and the non-centered
% cross-sectional covariate.
% =========================================================================
% FORMAT: swe_splitCovariate(cov, subject)
% -------------------------------------------------------------------------
% Inputs:
%  - cov: the covariate to split
%  - subject: vector contening the subject numbering 
% -------------------------------------------------------------------------
% Outputs: 
%  - crossCov: the centered cross-sectional covariate
%  - longCov: the longitudinal covariate
%  - avgCrossCov: the average of the time-varying covariate
% =========================================================================
% written by Bryan Guillaume
function [crossCov, longCov, avgCov] = swe_splitCovariate(cov, subject)

    uSubject = unique(subject);
    nSubj = length(uSubject);

    crossCov = 0 * cov;
    for i = 1:nSubj
      crossCov(subject == uSubject(i)) = mean(cov(subject == uSubject(i)));
    end
    longCov = cov - crossCov;
    avgCov = mean(crossCov);
    crossCov = crossCov - avgCov;

end
