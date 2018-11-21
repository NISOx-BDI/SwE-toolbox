function [crossCov, longCov, avgCov] = swe_splitCovariate(cov, subject)
% Splits a covariate into pure cross-sectional and longitudinal components.
% =========================================================================
% The split is done first by extracting the within-subect average to build 
% the cross-sectional covariate, which is also centered. Second, the 
% longitudinal covariate is simply the difference between the original 
% covariate and the non-centered cross-sectional covariate.
% =========================================================================
% FORMAT: [crossCov, longCov, avgCov] = swe_splitCovariate(cov, subject)
% -------------------------------------------------------------------------
% Inputs:
%  - cov: the covariate to split
%  - subject: vector contening the subject numbering 
% -------------------------------------------------------------------------
% Outputs: 
%  - crossCov: the centered cross-sectional covariate
%  - longCov: the longitudinal covariate
%  - avgCrossCov: the average value (scalar) of the time-varying covariate
% =========================================================================
% written by Bryan Guillaume
% Version Info:  $Format:%ci$ $Format:%h$

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
