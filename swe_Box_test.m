function [p,F,f1,f2] = swe_Box_test()
% Box's test of Compound Symmetry for the SwE toolbox
% FORMAT [p,F,f1,f2] = swe_Box_test()
%
% p     - vector of p-values
% F     - vector of F-scores
% f1    - numerator degrees of freedom
% f2    - denominator degrees of freedom
% 
%___________________________________________________________________________
%
% Note that the Box's test of Compound Symmetry assumes that there is no 
% missingness in the data. See Box (1950) for more information.

% Written by Bryan Guillaume
% $Id$

%select a SwE.
P     = cellstr(spm_select(1,'^SwE\.mat$','Select SwE.mat[s]'));
swd     = fileparts(P{1});
load(fullfile(swd,'SwE.mat'));

if SwE.Gr.nGr~=1
     error('More than 1 common covariance matrix in the design. Please use a design with only one common covariance matrix');
end

for i = 1:SwE.Subj.nSubj
    if length(SwE.Vis.iVis(SwE.Subj.iSubj==SwE.Subj.uSubj(i)))~= SwE.Vis.nVis
        error('Presence of missing data. Please provide a design without missing data');
    end
end

k       = SwE.Vis.nVis;
df      = SwE.dof.edof_Gr;
f1      = (k^2+k-4)/2;
A1      = k*(k+1)^2*(2*k-3)/(6*df*(k-1)*(k^2+k-4));
A2      = (k-1)*k*(k+1)*(k+2)/(6*df^2*(k^2+k-4));
f2      = (f1+2)/(A2-A1^2);
coeff   = df * (1-A1-f1/f2)/f1;
cov_vis = spm_get_data(SwE.Vcov_vis,SwE.xVol.XYZ);
p       = zeros(1,size(cov_vis,2));
F       = zeros(1,size(cov_vis,2));

str = 'CS SwE computation';
spm_progress_bar('Init',100,str,'');
for iVox = 1:size(cov_vis,2)
    S = zeros(k,k);
    it=0;
    for j = 1:k
        for jj = j:k
            it = it+1;
            S(j,jj) = cov_vis(it,iVox);
        end
    end
    S = S+S'-diag(diag(S));
    detS  = det(S);   
    var   = trace(S)/k;
    cov   = sum(sum(triu(S,1)))/k/(k-1)*2;
    detS0 = (var+(k-1)*cov)*(var-cov)^(k-1);
    F(iVox) = coeff*log(detS0/detS);% follow a F(f1,f2) under the Null: S is CS
    p(iVox) = -log10(fcdf(F(iVox),f1,f2,'upper'));
    spm_progress_bar('Set',100*(iVox/size(cov_vis,2)));
end
spm_progress_bar('Clear')

CS_test_F = struct(...
    'fname',  sprintf('swe_vox_Fstat-Box%s', spm_file_ext),...
    'dim',    SwE.xVol.DIM',...
    'dt',     [16, spm_platform('bigend')],...
    'mat',    SwE.xVol.M,...
    'pinfo',  [1,0,0]',...
    'descrip',sprintf('CS_test F-scores'));

CS_test_p = struct(...
    'fname',  sprintf('swe_vox_Fstat-Box_lp%s', spm_file_ext),...
    'dim',    SwE.xVol.DIM',...
    'dt',     [16, spm_platform('bigend')],...
    'mat',    SwE.xVol.M,...
    'pinfo',  [1,0,0]',...
    'descrip',sprintf('CS_test p-values'));

tmp           = zeros(SwE.xVol.DIM');
Q             = cumprod([1,SwE.xVol.DIM(1:2)'])*SwE.xVol.XYZ - ...
    sum(cumprod(SwE.xVol.DIM(1:2)'));
tmp(Q)        = F;
spm_write_vol(CS_test_F,tmp);
tmp(Q)        = p;
spm_write_vol(CS_test_p,tmp);



