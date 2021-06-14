function [groups,subjects,visits,X] = swe_GenFactorialDes(nVis,nSub)
% Create key input variables for a full factorial design
% =========================================================================
% Based on the number of visits and number of subjects per group, this
% function creates the four essential inputs for a SwE analysis:
% groups, subjects and visit vectors, and the design matrix X.
% 
% Crucially, the function assumes that visit indices vary fastest, and 
% groups vary slowest (i.e. the first scans are all the visits for first 
% subject in group 1.)
%
% =========================================================================
% FORMAT: [groups,subjects,visits,X] = swe_GenFactorialDes(nVis,nSub)
% -------------------------------------------------------------------------
% Inputs:
%  - nVis: Number of visits (or repeated measurements) per subject
%  - nSub: Number of subjects per group, a vector, one element per group
% -------------------------------------------------------------------------
% Outputs:
%  - groups: Groups indicator, values ranging from 1 to length(nSub)
%  - subjects: Subjects indicator, values ranging from 1 to sum(nSub)
%  - visits: Visits indicator, values ranging from 1 to nVis
%  - X: Design matrix, nVis*sum(nSum) rows, nvis*length(nSub) columns
% =========================================================================
% 
% Unbalanced designs can be obtained by carefully removing elements of 
% the groups, subjects & visits vectors and the corresponding rows of X.
%
% =========================================================================
% written by Thomas Nichols
% Version Info:  $Format:%ci$ $Format:%h$

[groups,visits,subjects,X]=deal([]);

Xv = eye(nVis);
for ig = 1:length(nSub)
    groups   = [groups,   repelem(ig,nVis*nSub(ig))];
    subjects = [subjects, repelem(1:nSub(ig),nVis)+sum(nSub(1:ig-1))];
    visits   = [visits,   repmat(1:nVis,[1 nSub(ig)])];

    X0 = kron(ones(nSub(ig),1),Xv);
    X  = [X                           zeros(size(X,1),size(X0,2));
          zeros(size(X0,1),size(X,2)) X0];
end

