function swe_cp_WB(SwE)

%-Say hello
%--------------------------------------------------------------------------
Finter = spm('CreateIntWin','off');
set(Finter,'name','SwE estimation');
set(Finter,'vis','on')

%-Change to SwE.swd if specified
%--------------------------------------------------------------------------
disp('0')
disp(SwE)
try
  cd(SwE.swd);
catch %#ok<*CTCH>
  SwE.swd = pwd;
end

disp('1')

%-Ensure data are assigned
%--------------------------------------------------------------------------
try
  SwE.xY.VY;
catch
  spm('alert!','Please assign data to this design', mfilename);
  spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
  return
end

disp('2')

%-Check if we have data in a.mat format and set some variables accordingly
%--------------------------------------------------------------------------
[~,~,file_ext] = fileparts(SwE.xY.P{1});
isMat          = strcmpi(file_ext,'.mat');

if ~isMat
    isMeshData = spm_mesh_detect(SwE.xY.VY);
    if isMeshData
        file_ext = '.gii';
    else
        file_ext = spm_file_ext;
    end
else
    isMeshData = false;
end

%-Delete files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(SwE.swd,sprintf('swe_vox_mask%s',file_ext)),'file') == 2
  
  str = {'Current directory contains SwE estimation files:',...
    'pwd = ',SwE.swd,...
    'Existing results will be overwritten!'};
  if spm_input(str,1,'bd','stop|continue',[1,0],1)
    spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
    return
  else
    warning('Overwriting old results\n\t (pwd = %s) ',SwE.swd); %#ok<WNTAG>
    try SwE = rmfield(SwE,'xVol'); end %#ok<TRYNC>
  end
end

disp('3')

files = {'^mask\..{3}$','^ResMS\..{3}$','^RPV\..{3}$',...
    '^beta_.{4}\..{3}$','^con_.{4}\..{3}$','^ResI_.{4}\..{3}$',...
    '^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$',...
    '^cov_beta_.{4}_.{4}\..{3}$', '^cov_vis_.{4}_.{4}_.{4}\..{3}$',...
    '^edf_.{4}\..{3}$', '^spm\w{4}_.{4}\..{3}$'};

for i = 1:length(files)
  j = spm_select('List',SwE.swd,files{i});
  for k = 1:size(j,1)
    spm_unlink(deblank(j(k,:)));
  end
end

%==========================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%==========================================================================

%-Initialise
%==========================================================================
fprintf('%-40s: %30s','Initialising parameters','...computing');        %-#
xX            = SwE.xX;
[nScan, nBeta] = size(xX.X);
nCov_beta     = (nBeta+1)*nBeta/2;
pX            = pinv(xX.X); % pseudo-inverse
iSubj         = SwE.Subj.iSubj;
uSubj         = unique(iSubj);
nSubj         = length(uSubj);

SwE.Subj.uSubj = uSubj;
SwE.Subj.nSubj = nSubj;

%-WB variables
%
WB = SwE.WB;
conWB = WB.con;
nSizeCon = size(conWB,1);
rankCon = rank(conWB);

% If clusterWise inference, force the U-SwE to be used
if WB.clusterWise == 1 && WB.RSwE ==1
  WB.RSwE = 0;
  SwE.WB.RSwE = 0;
  if spm_matlab_version_chk('7') >=0
    save('SwE','SwE','-V6');
  else
    save('SwE','SwE');
  end
end

% If clusterWise inference and .mat format, check for the presence of
% spatial information
if isMat && WB.clusterWise == 1
  if isfield(SwE.WB.clusterInfo, 'Vxyz')
    XYZ = importdata(SwE.WB.clusterInfo.Vxyz{1});
    if size(XYZ,1) ~=3 && size(XYZ,2) ~=3
      error('voxel coodinates do not seem correct')
    elseif size(XYZ,2) ==3
      XYZ = XYZ';
    end
  elseif isfield(SwE.WB.clusterInfo, 'Vfaces')
    faces = importdata(SwE.WB.clusterInfo.Vfaces{1});
    if size(faces,1) ~=3 && size(faces,2) ~=3
      error('faces coodinates do not seem correct')
    elseif size(faces,1) ==3
      faces = faces';
    end
  else
    error('clusterWise inference cannot be done without spatial information when inputs are in ".mat" format. Please supply faces coordinates (faces or tris) for surface data or voxel coordinates (XYZ_vox) for volumetric data');
  end
end

% small sample correction (for WB)
[corrWB, tmpR2] = swe_resid_corr(SwE, WB.RWB, WB.SS, pX);

% small sample correction (for parametric)
[corr, tmpR2] = swe_resid_corr(SwE, WB.RSwE, SwE.SS, pX, tmpR2);

%-detect if the design matrix is separable (a little bit messy, but seems to do the job)
%
iGr_dof   = zeros(1,nScan);
iBeta_dof = zeros(1,nBeta);
it = 0;
while ~all(iGr_dof)
    it = it + 1;
    scan = find(iGr_dof==0,1);
    if any(xX.X(scan,:)) % handle the case where a row is all 0s (BG - 05/08/2016; Thanks to Ged Ridgway for finding the bug)
      for i = find(iGr_dof==0)
          if any((xX.X(i,:) & xX.X(scan,:)))
              iGr_dof(i) = it;
          end
      end
    else
      iGr_dof(scan) = it;
    end
end
%need to check if the partition is correct
while 1
  uGr_dof = unique(iGr_dof);
  nGr_dof = length(uGr_dof);
  tmp = zeros(nGr_dof,nBeta);
  for i = 1:nGr_dof
    tmp(i,:) = any(xX.X(iGr_dof==i,:));
  end
  if nGr_dof==1 | all(sum(tmp, 1)==1) %#ok<OR2>
    break % all is ok, just stop the while
  else
    ind1 = find(sum(tmp, 1)>1,1); % detect the first column in common
    ind2 = find(tmp(:,ind1)==1); % detect the groups to be fused
    for ii = ind2'
      iGr_dof(iGr_dof==ii) = ind2(1); % fuse the groups
    end
  end
end
nSubj_dof = zeros(1,nGr_dof);
for i = 1:nGr_dof % renumber to avoid gaps in the numbering
  iGr_dof(iGr_dof==uGr_dof(i)) = i;
  iBeta_dof(tmp(i,:)==1) = i;
  nSubj_dof(i) = length(unique(iSubj(iGr_dof==uGr_dof(i))));
end

pB_dof   = zeros(1,nGr_dof); 
for i=1:nBeta    
    tmp=1;
    for ii=1:nSubj
        if length(unique(xX.X(iSubj==uSubj(ii)&iGr_dof'==iBeta_dof(i),i)))>1
            tmp=0;
            break
        end
    end
    if tmp == 1
        pB_dof(iBeta_dof(i)) = pB_dof(iBeta_dof(i)) + 1;
    end
end

disp('7')

%-effective dof for each subject
edof_Subj = zeros(1,nSubj);
for i = 1:nSubj
    edof_Subj(i) = 1 - pB_dof(iGr_dof(iSubj==uSubj(i)))/...
        nSubj_dof(iGr_dof(iSubj==uSubj(i)));
end

%-degrees of freedom estimation type
if isfield(SwE.type,'modified')
    dof_type = SwE.type.modified.dof_mo;
else
    dof_type = SwE.type.classic.dof_cl;        
end

if dof_type == 0 % so naive estimation is used
  if nSizeCon==1
    ind = find(conWB ~= 0);
  else
    ind = find(any(conWB~=0));
  end
  edf = sum(nSubj_dof(unique(iBeta_dof(ind))) - pB_dof(unique(iBeta_dof(ind))));
  
  dof_cov = zeros(1,nBeta);
  for i = 1:nBeta
      dof_cov(i) = nSubj_dof(iBeta_dof(i)) - ...
          pB_dof(iBeta_dof(i));    
  end
else
    edf = NaN;
end

%-preprocessing for the modified SwE
if isfield(SwE.type,'modified')
  iVis      = SwE.Vis.iVis;
  iGr       = SwE.Gr.iGr;
  uGr       = unique(iGr);
  nGr       = length(uGr);
  SwE.Gr.uGr       = uGr;
  SwE.Gr.nGr       = nGr;
  
  % info specific for each group
  uVis_g = cell(1,nGr); % unique visits for each group
  nVis_g = zeros(1,nGr); % number of visits for each group
  uSubj_g = cell(1,nGr); % unique subjects for each group
  nSubj_g = zeros(1,nGr); % number of subjects for each group
  for g = 1:nGr
    uVis_g{g}  = unique(iVis(iGr==uGr(g)));
    nVis_g(g)  = length(uVis_g{g});
    uSubj_g{g} = unique(iSubj(iGr==uGr(g)));
    nSubj_g(g) = length(uSubj_g{g});
  end
  nCov_vis_g  = nVis_g.*(nVis_g+1)/2; % number of covariance elements to be estimated for each group
  nCov_vis    = sum(nCov_vis_g); % total number of covariance elements to be estimated
  
  % Save nVis_g and uVis_g.
  SwE.Vis.uVis_g = uVis_g;
  SwE.Vis.nVis_g = nVis_g;
  
  % Flags matrices indicating which residuals have to be used for each covariance element
  Flagk  = false(nCov_vis,nScan); % Flag indicating scans corresponding to visit k for each covariance element
  Flagkk = false(nCov_vis,nScan); % Flag indicating scans corresponding to visit kk for each covariance element
  Ind_Cov_vis_diag     = nan(1,sum(nVis_g)); % index of the diagonal elements
  Ind_Cov_vis_off_diag = nan(1,nCov_vis - sum(nVis_g)); % index of the off-diagonal elements
  Ind_corr_diag=nan(nCov_vis,2); % index of the 2 corresponding diagonal elements
  iGr_Cov_vis_g = nan(1,nCov_vis);
  it = 0; it2 = 0; it3 = 0;
  for g = 1:nGr
    for k = 1:nVis_g(g)
      for kk = k:nVis_g(g)
        it = it + 1;
        id = intersect(iSubj(iGr==uGr(g) & iVis==uVis_g{g}(k)),...
          iSubj(iGr==uGr(g) & iVis==uVis_g{g}(kk))); % identifiaction of the subjects with both visits k & kk
        Flagk(it,:)  = ismember(iSubj,id) & iVis==uVis_g{g}(k);
        Flagkk(it,:) = ismember(iSubj,id) & iVis==uVis_g{g}(kk);
        if k==kk
          it2 = it2+1;
          it4 = it2;
          Ind_Cov_vis_diag(it2)     = it;
        else
          it3 = it3 + 1;
          it4 = it4 + 1;
          Ind_Cov_vis_off_diag(it3) = it;
        end
        Ind_corr_diag(it,:) = [it2 it4];
        iGr_Cov_vis_g(it) = g;
      end
    end
  end
  
  % Record igr_Cov_vis_g.
  SwE.WB.iGr_Cov_vis_g = iGr_Cov_vis_g;
  
  % weights for the vectorised SwE
  weight = NaN(nCov_beta,nCov_vis);

  it=0;
  for j = 1:nBeta
    for jj = j:nBeta
      it=it+1;
      for jjj = Ind_Cov_vis_diag
        weight(it,jjj) = pX(j,Flagk(jjj,:))*pX(jj,Flagk(jjj,:))';
      end
      for jjj = Ind_Cov_vis_off_diag
        weight(it,jjj) = pX(j,Flagk(jjj,:))*pX(jj,Flagkk(jjj,:))' + ...
          pX(j,Flagkk(jjj,:))*pX(jj,Flagk(jjj,:))';
      end
    end
  end
  % Weight giving only the contrasted SwE (WB)
  weightR = pinv(swe_duplication_matrix(nSizeCon)) * kron(conWB,conWB) * swe_duplication_matrix(nBeta) * weight; % used to compute the R SwE R' 
  Wg = cell(nGr,1);
  Wg_testII = cell(nGr,1);
  Wg_testIII = cell(nGr,1);
  tmp = eye(nSizeCon^2);  

  for g = 1:nGr
    Wg{g} = kron(weightR(:,iGr_Cov_vis_g==g),weightR(:,iGr_Cov_vis_g==g)) * swe_duplication_matrix(nCov_vis_g(g));
    Wg_testII{g} = sum(kron(swe_duplication_matrix(nSizeCon),swe_duplication_matrix(nSizeCon)), 1) * Wg{g};
    Wg_testIII{g} = tmp(:)' * (kron(swe_duplication_matrix(nSizeCon),swe_duplication_matrix(nSizeCon))) * Wg{g};
  end
   
  SwE.WB.Wg{1} = Wg;
  SwE.WB.Wg{2} = Wg_testII;
  SwE.WB.Wg{3} = Wg_testIII;
  
%-compute the effective dof from each homogeneous group if dof_type
    switch dof_type
      case 1
        dofMat = NaN;
        edof_Gr = zeros(1,nGr);
        nSubj_g = zeros(1,nGr);
        for g = 1:nGr
          nSubj_g(g) = length(unique(iSubj(iGr == g)));
          tmp = 0;
          for j = 1:nSubj_g(g)
            tmp = tmp + 1/edof_Subj(uSubj == uSubj_g{g}(j));
          end
          edof_Gr(g) = nSubj_g(g)^2/tmp;
        end
      case {2,3} % compute a matrix containing the variables linked to the degrees of freedom (for test II and III)
        dofMat = cell(nGr,1);
        for g = 1:nGr
          dofMat{g} = zeros(nCov_vis_g(g));
          it1 =0;
          for i  = 1:nVis_g(g)
            for j  = i:nVis_g(g)
              it1 = it1 + 1;
              it2 = 0;
              for a = 1:nVis_g(g)
                for b = a:nVis_g(g)
                  it2 = it2 + 1;
                  mij = 0;mab = 0;tmp = 0;
                  for ii = 1:nSubj_g(g)
                    mij = mij + 1*(...
                      any(iSubj==uSubj_g{g}(ii) & iVis==uVis_g{g}(i)) &...
                      any(iSubj==uSubj_g{g}(ii) & iVis==uVis_g{g}(j)));
                    mab = mab + 1*(...
                      any(iSubj==uSubj_g{g}(ii) & iVis==uVis_g{g}(a)) &...
                      any(iSubj==uSubj_g{g}(ii) & iVis==uVis_g{g}(b)));
                    tmp = tmp + 1*(...
                      any(iSubj==uSubj_g{g}(ii) & iVis==uVis_g{g}(a)) &...
                      any(iSubj==uSubj_g{g}(ii) & iVis==uVis_g{g}(b)) &...
                      any(iSubj==uSubj_g{g}(ii) & iVis==uVis_g{g}(i)) &...
                      any(iSubj==uSubj_g{g}(ii) & iVis==uVis_g{g}(j)))...
                      /edof_Subj(uSubj==uSubj_g{g}(ii));
                  end
                  dofMat{g}(it1,it2) = tmp/mij/mab;
                end
              end
            end
          end
          dofMat{g}(isnan(dofMat{g})) = 0;
        end
        clear tmp mij mab
    end
end

disp('10')

%-preprocessing for the classic SwE
if isfield(SwE.type,'classic')
  nVis_i        = zeros(1,nSubj);
  for i = 1:nSubj
    nVis_i(i) = sum(uSubj(i)==iSubj);
  end
  nCov_vis      = sum(nVis_i.*(nVis_i+1)/2); % total number of covariance elements to be estimated
  weight        = NaN(nCov_beta,nCov_vis);
  Ind_Cov_vis_classic = NaN(1,nCov_vis);
  Indexk  = NaN(1,nCov_vis);
  Indexkk = NaN(1,nCov_vis);
  
  it=0;
  for j = 1:nBeta
    for jj = j:nBeta
      it = it + 1;
      it2 = 0;
      for i = 1:nSubj
        ind_i=find(iSubj == uSubj(i));
        for ii = 1:nVis_i(i)
          it2 = it2 + 1;
          weight(it,it2) = pX(j,ind_i(ii))*pX(jj,ind_i(ii));
          Ind_Cov_vis_classic(it2) = i;
          Indexk(it2)  = ind_i(ii);
          Indexkk(it2) = ind_i(ii);
          for iii = (ii+1):nVis_i(i)
            it2 = it2 + 1;
            weight(it,it2) = pX(j,ind_i([ii,iii]))*pX(jj,ind_i([iii,ii]))';
            Ind_Cov_vis_classic(it2) = i;
            Indexk(it2)  = ind_i(ii);
            Indexkk(it2) = ind_i(iii);
          end
        end
      end
    end
  end
  %-compute the effective dof from each homogeneous group (here, subject)
%   if dof_type == 1
%     edof_Gr = edof_Subj;
%   end
  weightR = pinv(swe_duplication_matrix(nSizeCon)) * kron(conWB,conWB) * swe_duplication_matrix(nBeta) * weight; % used to compute the R SwE R' 
end

disp('11')

%-If xM is not a structure then assume it's a vector of thresholds
%--------------------------------------------------------------------------
try
  xM = SwE.xM;
catch
  xM = -Inf(nScan,1);
end
if ~isstruct(xM)
  xM = struct('T',    [],...
    'TH',   xM,...
    'I',    0,...
    'VM',   {[]},...
    'xs',   struct('Masking','analysis threshold'));
end

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');               %-#

if ~isMat
  %-Image dimensions and data
  %==========================================================================
  VY       = SwE.xY.VY;
  spm_check_orientations(VY);
  
  % check files exists and try pwd
  %--------------------------------------------------------------------------
  for i = 1:numel(VY)
    if ~spm_existfile(VY(i).fname)
      [p,n,e]     = fileparts(VY(i).fname);
      VY(i).fname = [n,e];
    end
  end
  
  M        = VY(1).mat;
  DIM      = VY(1).dim(1:3)';
  VOX      = sqrt(diag(M(1:3, 1:3)'*M(1:3, 1:3)))';
  xdim     = DIM(1); ydim = DIM(2); zdim = DIM(3);
  %vFWHM    = SwE.vFWHM; to be added later (for the variance smoothing)
  % YNaNrep  = spm_type(VY(1).dt(1),'nanrep');
  YNaNrep = VY(1).dt(2);
    
  fprintf('%-40s: %30s','Output images','...initialising');           %-#
  
  %-Initialise new mask name: current mask & conditions on voxels
  %----------------------------------------------------------------------
  VM    = struct('fname',  sprintf('swe_vox_mask%s', file_ext),...
    'dim',    DIM',...
    'dt',     [spm_type('uint8') spm_platform('bigend')],...
    'mat',    M,...
    'pinfo',  [1 0 0]',...
    'descrip','swe_cp_WB:resultant analysis mask');
  VM    = spm_create_vol(VM);
  
  %-Initialise original parametric score image, T or F
  %----------------------------------------------------------------------
  if WB.stat=='T'
    eSTAT='z';
  else % F stat
    eSTAT='x';
  end
  
  Vscore = swe_create_vol(sprintf('swe_vox_%cstat_c%02d%s', WB.stat, 1, file_ext), DIM, M,...
			  sprintf('Original parametric %c statistic data.', WB.stat));
  
  %-Initialise parametric P-Value image
  %----------------------------------------------------------------------
  
  VlP = swe_create_vol(sprintf('swe_vox_%cstat_lp_c%02d%s', WB.stat, 1, file_ext), DIM, M,...
                      'Original parametric -log10(P) value data (positive).');
  
  if WB.stat=='T'
        VlP_Neg = swe_create_vol(sprintf('swe_vox_%cstat_lp_c%02d%s', WB.stat, 2, file_ext), DIM, M,...
                               'Original parametric -log10(P) value data (negative).');
  end
  
  %-Initialise converted parametric score image
  %----------------------------------------------------------------------
  VcScore = swe_create_vol(sprintf('swe_vox_%c%cstat_c%02d%s', eSTAT, WB.stat, 1, file_ext), DIM, M,...
			   sprintf('Parametric %c statistic data derived from %c-Statistic data.', eSTAT, WB.stat));
  
  %-Initialise residual images for the resampling
  %----------------------------------------------------------------------
  
  for i = 1:nScan
      if WB.RWB == 1
        descrip = sprintf('adjusted restricted residuals (%04d)', i);
      else
        descrip = sprintf('adjusted unrestricted residuals (%04d)', i);
      end
      VResWB(i) = swe_create_vol(sprintf('swe_vox_resid_y%04d%s', i, file_ext), DIM, M, descrip);
  end
  
  fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');    %-#
  
  %-Initialise fitted data images for the resampling
  %----------------------------------------------------------------------
  
  for i = 1:nScan
      if WB.RWB == 1
         descrip = sprintf('restricted fitted data  (%04d)', i);
      else
         descrip = sprintf('unrestricted fitted data (%04d)', i);
      end
      VYWB(i) = swe_create_vol(sprintf('swe_vox_fit_y%04d%s',i,file_ext), DIM, M, descrip);
  end
  
  %-Initialise result images
  %----------------------------------------------------------------------
  VlP_wb_pos = swe_create_vol(sprintf('swe_vox_%cstat_lp-WB_c%02d%s', WB.stat, 1, file_ext), DIM, M,...
                              '-log10(uncor. non-para. P, +ve)');

  VlP_wb_FWE_pos = swe_create_vol(sprintf('swe_vox_%cstat_lpFWE-WB_c%02d%s', WB.stat, 1, file_ext), DIM, M,...
                                  '-log10(FWE-corr. P, +ve)');
  
  VlP_wb_FDR_pos = swe_create_vol(sprintf('swe_vox_%cstat_lpFDR-WB_c%02d%s', WB.stat, 1, file_ext), DIM, M,...
                                  '-log10(FDR-corr. P, +ve)');

  if WB.stat=='T'
    VlP_wb_neg = swe_create_vol(sprintf('swe_vox_%cstat_lp-WB_c%02d%s', WB.stat, 2, file_ext), DIM, M,...
                              '-log10(uncor. non-para. P, -ve)');
    
    VlP_wb_FWE_neg = swe_create_vol(sprintf('swe_vox_%cstat_lpFWE-WB_c%02d%s', WB.stat, 2, file_ext), DIM, M,...
                                     '-log10(FWE-corr. P, -ve)');
    
    VlP_wb_FDR_neg = swe_create_vol(sprintf('swe_vox_%cstat_lpFDR-WB_c%02d%s', WB.stat, 2, file_ext), DIM, M,...
                                     '-log10(FDR-corr. P, -ve)');

  end
  
  % Converted score for WB.
  VcScore_wb_pos = swe_create_vol(sprintf('swe_vox_%c%cstat-WB_c%02d%s', eSTAT, WB.stat, 1, file_ext), DIM, M,...
                                 'Z score image for wild bootstrap voxelwise results.');
  
  if WB.clusterWise == 1
      
    % We also need cluster p value maps here.
    VlP_wb_clusterFWE_pos = swe_create_vol(sprintf('swe_clustere_%cstat_lpFWE-WB_c%02d%s', WB.stat, 1, file_ext), DIM, M,...
                                           '-log10(clusterFWE-corr. P, +ve)');
    
    if WB.stat=='T'
      VlP_wb_clusterFWE_neg = swe_create_vol(sprintf('swe_clustere_%cstat_lpFWE-WB_c%02d%s', WB.stat, 2, file_ext), DIM, M,...
                                             '-log10(clusterFWE-corr. P, -ve)');
    end
  end
  
  fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');    %-#
  %==========================================================================
  % - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
  %==========================================================================
  
  %-MAXMEM is the maximum amount of data processed at a time (bytes)
  %--------------------------------------------------------------------------
  MAXMEM = spm_get_defaults('stats.maxmem');
  mmv    = MAXMEM/8/nScan;
  blksz  = min(xdim*ydim,ceil(mmv));                             %-block size
  nbch   = ceil(xdim*ydim/blksz);                                %-# blocks
  nbz    = max(1,min(zdim,floor(mmv/(xdim*ydim)))); nbz = 1;     %-# planes forced to 1 so far
  blksz  = blksz * nbz;
  
  %-Initialise variables used in the loop
  %==========================================================================
  [xords, yords] = ndgrid(1:xdim, 1:ydim);
  xords = xords(:)'; yords = yords(:)';           % plane X,Y coordinates
  S     = 0;                                      % Volume (voxels)
  
  %-Initialise XYZ matrix of in-mask voxel co-ordinates (real space)
  %--------------------------------------------------------------------------
  XYZ   = zeros(3,xdim*ydim*zdim);
  
disp('13')
  %-Cycle over bunches blocks within planes to avoid memory problems
  %==========================================================================
  str   = 'parameter estimation';
  spm_progress_bar('Init',100,str,'');
  
  % activated voxels for cluster-wise inference
  if (WB.clusterWise == 1)
    activatedVoxels = false(0);
    maxClusterSize = nan(1, WB.nB + 1);
    activatedVoxelsNeg = false(0);
    if (WB.stat == 'T')
      maxClusterSizeNeg = nan(1, WB.nB + 1);
    end
  end
  maxScore = nan(1, WB.nB + 1);
  if (WB.stat == 'T')
    minScore = nan(1, WB.nB + 1);
  end
  
  for z = 1:nbz:zdim                       %-loop over planes (2D or 3D data)
    
    % current plane-specific parameters
    %----------------------------------------------------------------------
    CrPl         = z:min(z+nbz-1,zdim);       %-plane list
    zords        = CrPl(:)*ones(1,xdim*ydim); %-plane Z coordinates
    CrScore      = [];                        %-scores
    CrYWB         = [];                       %-fitted data under H0
    CrResWB       = [];                       %-residuals
    CrP          = [];                        %-parametric p-values
    if (WB.stat == 'T')
     CrPNeg       = [];                       %-negative parametric p-values
    end
    CrConScore   = [];                        %-converted score values. 
                                              % i.e. Z/X from T/F
    Q            = [];                        %-in mask indices for this plane
    
    for bch = 1:nbch                     %-loop over blocks
      
      %-Print progress information in command window
      %------------------------------------------------------------------
      if numel(CrPl) == 1
        str = sprintf('Plane %3d/%-3d, block %3d/%-3d',...
          z,zdim,bch,nbch);
      else
        str = sprintf('Planes %3d-%-3d/%-3d',z,CrPl(end),zdim);
      end
      if z == 1 && bch == 1
        str2 = '';
      else
        str2 = repmat(sprintf('\b'),1,72);
      end
      fprintf('%s%-40s: %30s',str2,str,' ');
      
      
      %-construct list of voxels in this block
      %------------------------------------------------------------------
      I     = (1:blksz) + (bch - 1)*blksz;       %-voxel indices
      I     = I(I <= numel(CrPl)*xdim*ydim);     %-truncate
      xyz   = [repmat(xords,1,numel(CrPl)); ...
        repmat(yords,1,numel(CrPl)); ...
        reshape(zords',1,[])];
      xyz   = xyz(:,I);                          %-voxel coordinates
      nVox  = size(xyz,2);                       %-number of voxels
      
      %-Get data & construct analysis mask
      %=================================================================
      fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...read & mask data')
      Cm    = true(1,nVox);                      %-current mask
      
      %-Compute explicit mask
      % (note that these may not have same orientations)
      %------------------------------------------------------------------
      for i = 1:length(xM.VM)
        
        %-Coordinates in mask image
        %--------------------------------------------------------------
        j = xM.VM(i).mat\M*[xyz;ones(1,nVox)];
        
        %-Load mask image within current mask & update mask
        %--------------------------------------------------------------
	disp('Before')
        Cm(Cm) = spm_get_data(xM.VM(i),j(:,Cm),false) > 0;
	disp('After')
      end
      
      %-Get the data in mask, compute threshold & implicit masks
      %------------------------------------------------------------------
      Y     = zeros(nScan,nVox);
      for i = 1:nScan
        
        %-Load data in mask
        %--------------------------------------------------------------
        if ~any(Cm), break, end                %-Break if empty mask
        Y(i,Cm)  = spm_get_data(VY(i),xyz(:,Cm),false);
        
        Cm(Cm)   = Y(i,Cm) > xM.TH(i);         %-Threshold (& NaN) mask
        if xM.I && ~YNaNrep && xM.TH(i) < 0    %-Use implicit mask
          Cm(Cm) = abs(Y(i,Cm)) > eps;
        end
      end
      
      %-Mask out voxels where data is constant in at least one separable
      % matrix design either in a visit category or within-subject (BG - 27/05/2016)
      %------------------------------------------------------------------
      [Cm, Y, CrS] = swe_mask_seperable(SwE, Cm, Y, iGr_dof);
      
      %==================================================================
      %-Proceed with General Linear Model (if there are voxels)
      %==================================================================
      if CrS
        
        %-General linear model: Ordinary least squares estimation
        %--------------------------------------------------------------
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...estimation');%-#
        
        beta  = pX*Y;                     %-Parameter estimates
        
        % restricted fitted data
        if WB.RWB == 1
            [resWB, YWB]=swe_fit(SwE, Y, tmpR2, corrWB, beta, SwE.WB.SS);
        else 
            [resWB, YWB]=swe_fit(SwE, Y, xX.X, corrWB, beta, SwE.WB.SS);
        end

        if WB.RSwE == 1
            res=swe_fit(SwE, Y, tmpR2, corr, beta, SwE.SS);
        else 
            res=swe_fit(SwE, Y, xX.X, corr, beta, SwE.SS);
        end

        clear Y                           %-Clear to save memory
        %-Estimation of the data variance-covariance components (modified SwE)
        %-SwE estimation (classic version)
        %--------------------------------------------------------------
        if isfield(SwE.type,'modified')
          Cov_vis=zeros(nCov_vis,CrS);
          for i = Ind_Cov_vis_diag
            Cov_vis(i,:) = mean(res(Flagk(i,:),:).^2, 1);
          end
          
          % Check if some voxels have variance < eps and mask them
          tmp = ~any(Cov_vis(Ind_Cov_vis_diag,:) < eps); % modified by BG on 29/08/16
          if any(~tmp)
            beta    = beta(:,tmp);
            resWB   = resWB(:,tmp);
            res     = res(:,tmp);
            YWB      = YWB(:,tmp);
            Cm(Cm)  = tmp;
            CrS     = sum(Cm);
            Cov_vis = Cov_vis(:,tmp);
          end
          if CrS % Check if there is at least one voxel left
            for i = Ind_Cov_vis_off_diag
              if any(Flagk(i,:))
                Cov_vis(i,:)= sum(res(Flagk(i,:),:).*res(Flagkk(i,:),:), 1).*...
                  sqrt(Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,1)),:).*...
                  Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,2)),:)./...
                  sum(res(Flagk(i,:),:).^2, 1)./...
                  sum(res(Flagkk(i,:),:).^2, 1));
              end
            end
            %NaN may be produced in cov. estimation when one correspondant
            %variance are = 0, so set them to 0
            Cov_vis(isnan(Cov_vis))=0;
            %need to check if the eigenvalues of Cov_vis matrices are >=0
            for g = 1:nGr
              for iVox = 1:CrS
                tmp = zeros(nVis_g(g));
                tmp(tril(ones(nVis_g(g)))==1) = Cov_vis(iGr_Cov_vis_g==g,iVox);
                tmp = tmp + tmp' - diag(diag(tmp));
                [V D] = eig(tmp);
                if any (diag(D)<0) %Bug corrected (BG - 19/09/13)
                  D(D<0) = 0;
                  tmp = V * D * V';
                  Cov_vis(iGr_Cov_vis_g==g,iVox) = tmp(tril(ones(nVis_g(g)))==1); %Bug corrected (BG - 19/09/13)
                end
              end
            end
          end
          cCovBc = weightR * Cov_vis;
        else
          cCovBc = 0;
          for i = 1:nSubj
            Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
              (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
            cCovBc = cCovBc + Cov_beta_i_tmp;
          end
        end
        
        % compute the score
        if (SwE.WB.stat == 'T')
          
          score = (conWB * beta) ./ sqrt(cCovBc);
          
          % hypothesis test, using clusterwise threshold if available.
          if (SwE.WB.clusterWise == 1)
            [p, activatedVoxels, activatedVoxelsNeg]=swe_hyptest(SwE, score, CrS, edf, cCovBc, Cov_vis, dofMat, activatedVoxels, activatedVoxelsNeg);
            clear CovcCovBc cCovBc
          else
            p=swe_hyptest(SwE, score, CrS, edf, cCovBc, Cov_vis, dofMat);
          end
          
          minScore(1) = min(minScore(1), min(score));
        else
          % need to loop at every voxel
          cBeta = conWB * beta;
          score = zeros(1, CrS);
          for iVox = 1:CrS
            cCovBc_vox = zeros(nSizeCon);
            cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
            cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
            score(iVox) = cBeta(:,iVox)' / cCovBc_vox * cBeta(:,iVox);
          end
          score = score / rankCon;
          
          % hypothesis test, using clusterwise threshold if available.
          if (SwE.WB.clusterWise == 1)
            [p, activatedVoxels]=swe_hyptest(SwE, score, CrS, edf, cCovBc, Cov_vis, dofMat, activatedVoxels);
          else
            p=swe_hyptest(SwE, score, CrS, edf, cCovBc, Cov_vis, dofMat);
          end
        end
        
        maxScore(1) = max(maxScore(1), max(score));
        
        %-Save betas etc. for current plane as we go along
        %----------------------------------------------------------
        CrYWB             = [CrYWB,    YWB]; %#ok<AGROW>
        CrResWB           = [CrResWB,  resWB]; %#ok<AGROW>
        CrScore           = [CrScore,  score]; %#ok<AGROW>
        CrP               = [CrP,      -log10(1-p)]; %#ok<AGROW>
        if (SwE.WB.stat == 'T')
            CrConScore    = [CrConScore, swe_invNcdf(p)]; %#ok<AGROW>
            CrPNeg        = [CrPNeg,   -log10(p)]; %#ok<AGROW>
        end
        if(SwE.WB.stat == 'F')
            CrConScore    = [CrConScore, spm_invXcdf(p, 1)]; %#ok<AGROW>
        end
        
      end % (CrS)
      
      %-Append new inmask voxel locations and volumes
      %------------------------------------------------------------------
      XYZ(:,S + (1:CrS)) = xyz(:,Cm);     %-InMask XYZ voxel coords
      Q                  = [Q I(Cm)];     %#ok<AGROW> %-InMask XYZ voxel indices
      S                  = S + CrS;       %-Volume analysed (voxels)
      
    end % (bch)
    
disp('15')
    %-Plane complete, write plane to image files (unless 1st pass)
    %======================================================================
    
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...saving plane'); %-#
    
    jj = NaN(xdim,ydim,numel(CrPl));
    
    %-Write Mask image
    %------------------------------------------------------------------
    if ~isempty(Q), jj(Q) = 1; end
	disp('lol')
    VM    = spm_write_plane(VM, ~isnan(jj), CrPl);
	disp('lol2')    

    %-Write WB fitted data images
    %------------------------------------------------------------------
    for i = 1:nScan
      if ~isempty(Q), jj(Q) = CrYWB(i,:); end
      VYWB(i) = spm_write_plane(VYWB(i), jj, CrPl);
    end
    
    %-Write WB residuals
    %------------------------------------------------------------------
    for i = 1:nScan
      if ~isempty(Q), jj(Q) = CrResWB(i,:); end
      VResWB(i) = spm_write_plane(VResWB(i), jj, CrPl);
    end
    
    %-Write parametric score image of the original data
    %------------------------------------------------------------------
    if ~isempty(Q), jj(Q) = CrScore; end
    Vscore = spm_write_plane(Vscore, jj, CrPl);
    
    %-Write parametric p-value image
    %------------------------------------------------------------------
    if ~isempty(Q), jj(Q) = CrP; end
    VlP = spm_write_plane(VlP, jj, CrPl);
    
    if WB.stat=='T'
        if ~isempty(Q), jj(Q) = CrPNeg; end
        VlP_Neg = spm_write_plane(VlP_Neg, jj, CrPl);
    end
    
    %-Write converted parametric score image of the original data
    %------------------------------------------------------------------
    if ~isempty(Q), jj(Q) = CrConScore; end
    VcScore = spm_write_plane(VcScore, jj, CrPl);
    
    %-Report progress
    %----------------------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done');
    spm_progress_bar('Set',100*(bch + nbch*(z - 1))/(nbch*zdim));
    
  end % (for z = 1:zdim)
  
  
  fprintf('\n');                                                          %-#
  spm_progress_bar('Clear')
  clear beta res Cov_vis CrScore CrYWB CrResWB Q jj%-Clear to save memory
  
  XYZ   = XYZ(:,1:S); % remove all the data not used
  
  % compute the max cluster size if needed (so many ways this can be
  % done... Not sure this solution is the best)
  if (SwE.WB.clusterWise == 1)
    LocActivatedVoxels = XYZ(:,activatedVoxels);
    clusterAssignment = spm_clusters(LocActivatedVoxels);
    nCluster     = max(clusterAssignment);
    clusterSize = histc(clusterAssignment,1:nCluster);
    
    if isempty(clusterSize)
      warning('no clusters survived the cluster-forming thresholding of the original data for positive effects!')
      maxClusterSize(1) = 0;
    else
      maxClusterSize(1) = max(clusterSize);
    end
    if (SwE.WB.stat == 'T')
      LocActivatedVoxelsNeg = XYZ(:,activatedVoxelsNeg);
      clusterAssignmentNeg = spm_clusters(LocActivatedVoxelsNeg);
      nClusterNeg     = max(clusterAssignmentNeg);
      clusterSizeNeg = histc(clusterAssignmentNeg,1:nClusterNeg);
      if isempty(clusterSizeNeg)
        warning('no clusters survived the cluster-forming thresholding of the original data for negative effects!')
        maxClusterSizeNeg(1) = 0;
      else
        maxClusterSizeNeg(1) = max(clusterSizeNeg);
      end
    end
  end
else % ".mat" format
  % check how the data image treat 0 (as NaN or not)
  YNaNrep = 0;
    
  fprintf('%-40s: %30s','Output images','...initialising');           %-#
  disp('16')
  
  fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');    %-#
  %==========================================================================
  % - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
  %==========================================================================
  
  
  %-Cycle over bunches blocks within planes to avoid memory problems
  %==========================================================================
  str   = 'parameter estimation';
  spm_progress_bar('Init',100,str,'');
  
  % activated voxels for cluster-wise inference
  if (WB.clusterWise == 1)
    activatedVoxels = false(0);
    maxClusterSize = nan(1, WB.nB + 1);
    if (WB.stat == 'T')
      activatedVoxelsNeg = false(0);
      maxClusterSizeNeg = nan(1, WB.nB + 1);
    end
  end
  maxScore = nan(1, WB.nB + 1);
  if (WB.stat == 'T')
    minScore = nan(1, WB.nB + 1);
  end
  
  %-Get data & construct analysis mask
  %=================================================================
  fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...read & mask data')
  
  %-Get the data in mask, compute threshold & implicit masks
  %------------------------------------------------------------------
  Y     = importdata(SwE.xY.P{1});
  % do some checking
  if ~isnumeric(Y)
    error('The input data is not a matrix. Please revised the model specification and specify a ".mat" file containing a unique matrix as data input file')
  elseif size(Y, 1) ~= SwE.nscan
    error('The input data does not have %i rows and thus is not compatible with the other specified variables. Please revised the model specification.', SwE.nscan)
  end
  
  nVox = size(Y, 2);
  
  %-Produce the mask
  Cm = true(1, nVox);
  %-Use the explicit mask if specified
  if length(SwE.xM.VM) == 1
    Cm(:) = importdata(SwE.xM.VM{1}) > 0;
  end
  %-check if some data need to be masked
  for i = 1:nScan
    if ~any(Cm), break, end                %-Break if empty mask
    Cm(Cm)   = Y(i,Cm) > xM.TH(i);         %-Threshold (& NaN) mask
    if xM.I && ~YNaNrep && xM.TH(i) < 0    %-Use implicit mask
      Cm(Cm) = abs(Y(i,Cm)) > eps;
    end
  end
  %-Mask out voxels where data is constant in at least one separable
  % matrix design either in a visit category or within-subject (BG - 27/05/2016)
  %------------------------------------------------------------------
  [Cm,Y,CrS] = swe_mask_seperable(SwE, Cm, Y, iGr_dof);
  
  if isfield(SwE.WB.clusterInfo, 'Vxyz')
    XYZ   = XYZ(:,Cm);
  end
  %==================================================================
  %-Proceed with General Linear Model (if there are voxels)
  %==================================================================
  if CrS
    
    %-General linear model: Ordinary least squares estimation
    %--------------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...estimation');%-#
    
    beta  = pX*Y;                     %-Parameter estimates
    
    if WB.RWB == 1
        [resWB, YWB]=swe_fit(SwE, Y, tmpR2, corrWB, beta, SwE.WB.SS);
    else 
        [resWB, YWB]=swe_fit(SwE, Y, xX.X, corrWB, beta, SwE.WB.SS);
    end
    
    if WB.RSwE == 1
        res=swe_fit(SwE, Y, tmpR2, corr, beta, SwE.SS);
    else 
        res=swe_fit(SwE, Y, xX.X, corr, beta, SwE.SS);
    end
    
    clear Y                           %-Clear to save memory
    %-Estimation of the data variance-covariance components (modified SwE)
    %-SwE estimation (classic version)
    %--------------------------------------------------------------
    if isfield(SwE.type,'modified')
      Cov_vis=zeros(nCov_vis,CrS);
      for i = Ind_Cov_vis_diag
        Cov_vis(i,:) = mean(res(Flagk(i,:),:).^2, 1);
      end
      
      % Check if some voxels have variance < eps and mask them
      tmp = ~any(Cov_vis(Ind_Cov_vis_diag,:) < eps); % modified by BG on 29/08/16
      if any(~tmp)
        beta    = beta(:,tmp);
        resWB   = resWB(:,tmp);
        res     = res(:,tmp);
        YWB     = YWB(:,tmp);
        Cm(Cm)  = tmp;
        CrS     = sum(Cm);
        Cov_vis = Cov_vis(:,tmp);
      end
      if CrS % Check if there is at least one voxel left
        for i = Ind_Cov_vis_off_diag
          if any(Flagk(i,:))
            Cov_vis(i,:)= sum(res(Flagk(i,:),:).*res(Flagkk(i,:),:), 1).*...
              sqrt(Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,1)),:).*...
              Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,2)),:)./...
              sum(res(Flagk(i,:),:).^2, 1)./...
              sum(res(Flagkk(i,:),:).^2, 1));
          end
        end
        %NaN may be produced in cov. estimation when one correspondant
        %variance are = 0, so set them to 0
        Cov_vis(isnan(Cov_vis))=0;
        %need to check if the eigenvalues of Cov_vis matrices are >=0
        for g = 1:nGr
          for iVox = 1:CrS
            tmp = zeros(nVis_g(g));
            tmp(tril(ones(nVis_g(g)))==1) = Cov_vis(iGr_Cov_vis_g==g,iVox);
            tmp = tmp + tmp' - diag(diag(tmp));
            [V, D] = eig(tmp);
            if any (diag(D)<0) %Bug corrected (BG - 19/09/13)
              D(D<0) = 0;
              tmp = V * D * V';
              Cov_vis(iGr_Cov_vis_g==g,iVox) = tmp(tril(ones(nVis_g(g)))==1); %Bug corrected (BG - 19/09/13)
            end
          end
        end
      end
      cCovBc = weightR * Cov_vis;
    else
      cCovBc = 0;
      for i = 1:nSubj
        Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
          (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
        cCovBc = cCovBc + Cov_beta_i_tmp;
      end
    end
    
    % compute the score
    if (SwE.WB.stat == 'T')
      
      score = (conWB * beta) ./ sqrt(cCovBc);
      
      if (SwE.WB.clusterWise == 1)
        [p, activatedVoxels, activatedVoxelsNeg]=swe_hyptest(SwE, score, CrS, edf, cCovBc, Cov_vis, dofMat, activatedVoxels, activatedVoxelsNeg);
        clear CovcCovBc cCovBc
      end
      
      minScore(1) = min(score);
    else
      % need to loop at every voxel
      cBeta = conWB * beta;
      score = zeros(1, CrS);
      for iVox = 1:CrS
        cCovBc_vox = zeros(nSizeCon);
        cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
        cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
        score(iVox) = cBeta(:,iVox)' / cCovBc_vox * cBeta(:,iVox);
      end
      score = score / rankCon;
      % Perform hypothesis test for activated regions.
      if (SwE.WB.clusterWise == 1)
        [p, activatedVoxels] = swe_hyptest(SwE, score, CrS, edf, cCovBc, Cov_vis, dofMat, activatedVoxels);
      end
    end
    maxScore(1) = max(score);
    
  end % (CrS)
  M           = [];
  DIM         = [];
  S           = CrS;  
  VM          = sprintf('swe_vox_mask%s', file_ext);
  Vscore      = sprintf('swe_vox_%cstat_c%02d%s', WB.stat, 1, file_ext);

  mask = Cm;       
  save(sprintf('swe_vox_mask%s',  file_ext), 'mask');
  clear mask
  
  tmp = score;
  score = nan(1, nVox);
  if (SwE.WB.stat == 'T')
      VT(:,Cm) = tmp;
      save(Vscore, 'VT');
      score = tmp;
      clear tmp VT
  else
      VF(:,Cm) = tmp;
      save(Vscore, 'VF');
      score = tmp;
      clear tmp VF
  end
  
  VlP = nan(1, nVox);
  VlP(:,Cm) = -log10(1-p);
  save(sprintf('swe_vox_%cstat_lp_c%02d%s', WB.stat, 1, file_ext), 'VlP');
  clear VlP

  if (SwE.WB.stat == 'T')
      
       VlP_neg = nan(1, nVox);
       VlP_neg(:,Cm) =  -log10(p);
       save(sprintf('swe_vox_%cstat_lp_c%02d%s', WB.stat, 2, file_ext), 'VlP_neg');
       clear VlP_neg
       
       z_map = nan(1, nVox);
       VZ(:,Cm) =  swe_invNcdf(p);
       save(sprintf('swe_vox_z%cstat_c%02d%s', WB.stat, 1, file_ext), 'VZ');
       clear VZ
  
  else
      
       x_map = nan(1, nVox);
       VX(:,Cm) =  spm_invXcdf(p, 1);
       save(sprintf('swe_vox_x%cstat_c%02d%s', WB.stat, 1, file_ext), 'VX');
       clear VX
       
  end
  
  fprintf('\n');                                                        %-#
  spm_progress_bar('Clear')
  clear res Cov_vis jj%-Clear to save memory
    
  % compute the max cluster size if needed (so many ways this can be
  % done... Not sure this solution is the best)
  if (SwE.WB.clusterWise == 1)
    if isfield(SwE.WB.clusterInfo, 'Vxyz')
      LocActivatedVoxels = XYZ(:,activatedVoxels);
      clusterAssignment = spm_clusters(LocActivatedVoxels);
    else %surface data      
      LocActivatedVoxels = false(nVox,1);
      LocActivatedVoxels(Cm) = activatedVoxels;
      clusterAssignment = spm_mesh_clusters(faces,LocActivatedVoxels);
      clusterAssignment = clusterAssignment(LocActivatedVoxels)';
      if isnan(clusterAssignment)
        clusterAssignment = [];
      end
    end
    nCluster     = max(clusterAssignment);
    clusterSize = histc(clusterAssignment,1:nCluster);
    
    if isempty(clusterSize)
      warning('no clusters survived the cluster-forming thresholding of the original data for positive effects!')
      maxClusterSize(1) = 0;
    else
      maxClusterSize(1) = max(clusterSize);
    end
    if (SwE.WB.stat == 'T')
      if isfield(SwE.WB.clusterInfo, 'Vxyz')
        LocActivatedVoxelsNeg = XYZ(:,activatedVoxelsNeg);
        clusterAssignmentNeg = spm_clusters(LocActivatedVoxelsNeg);
      else %surface data
        LocActivatedVoxelsNeg = false(nVox,1);
        LocActivatedVoxelsNeg(Cm) = activatedVoxelsNeg;
        clusterAssignmentNeg = spm_mesh_clusters(faces,LocActivatedVoxelsNeg);
        clusterAssignmentNeg = clusterAssignmentNeg(LocActivatedVoxelsNeg)';
        if isnan(clusterAssignmentNeg)
          clusterAssignmentNeg = [];
        end
      end
      nClusterNeg     = max(clusterAssignmentNeg);
      clusterSizeNeg = histc(clusterAssignmentNeg,1:nClusterNeg);
      if isempty(clusterSizeNeg)
        warning('no clusters survived the cluster-forming thresholding of the original data for negative effects!')
        maxClusterSizeNeg(1) = 0;
      else
        maxClusterSizeNeg(1) = max(clusterSizeNeg);
      end
    end
  end  
end 
%==========================================================================
% - P O S T   E S T I M A T I O N   C L E A N U P
%==========================================================================
if S == 0, spm('alert!','No inmask voxels - empty analysis!'); return; end

%-Save remaining results files and analysis parameters
%==========================================================================
fprintf('%-40s: %30s','Saving results','...writing');

%-place fields in SwE
%--------------------------------------------------------------------------
if isfield(SwE.WB, 'clusterInfo') && isfield(SwE.WB.clusterInfo, 'Vfaces')
  XYZ = [];
end

SwE.xVol.XYZ   = XYZ;               %-InMask XYZ coords (voxels)
SwE.xVol.M     = M;                 %-voxels -> mm
SwE.xVol.iM    = inv(M);            %-mm -> voxels
SwE.xVol.DIM   = DIM;               %-image dimensions
SwE.xVol.S     = S;
SwE.xVol.units = {'mm' 'mm' 'mm'};

if ~isMat
  SwE.WB.VYWB       = VYWB;               %-Filehandle - fitted data under H0
  SwE.WB.VResWB     = VResWB;             %-Filehandle - adjusted resticted residuals
  SwE.WB.Vscore     = Vscore;            %-Filehandle - score original image
end
SwE.WB.weightR    = weightR;
SwE.WB.corrWB     = corrWB;
SwE.WB.corr       = corr;

if SwE.WB.RWB == 1 || SwE.WB.RSwE == 1
  SwE.WB.tmpR2      = tmpR2;
end

% cluster-wise specific fields if needed
if (SwE.WB.clusterWise == 1)
  SwE.WB.clusterInfo.LocActivatedVoxels = LocActivatedVoxels;
  SwE.WB.clusterInfo.nCluster = nCluster;
  SwE.WB.clusterInfo.clusterAssignment = clusterAssignment;
  SwE.WB.clusterInfo.maxClusterSize = maxClusterSize;
  SwE.WB.clusterInfo.clusterSize = clusterSize;
  if (SwE.WB.stat == 'T')
    SwE.WB.clusterInfo.LocActivatedVoxelsNeg = LocActivatedVoxelsNeg;
    SwE.WB.clusterInfo.nClusterNeg = nClusterNeg;
    SwE.WB.clusterInfo.clusterAssignmentNeg = clusterAssignmentNeg;
    SwE.WB.clusterInfo.maxClusterSizeNeg = maxClusterSizeNeg;
    SwE.WB.clusterInfo.clusterSizeNeg = clusterSizeNeg;  
  end
end

SwE.WB.pX         = pX;
if isfield(SwE.type,'modified')
  SwE.WB.Ind_Cov_vis_diag = Ind_Cov_vis_diag;
  SwE.WB.Ind_Cov_vis_off_diag = Ind_Cov_vis_off_diag;
  SwE.WB.Ind_corr_diag = Ind_corr_diag;
  SwE.WB.Flagk      = Flagk;
  SwE.WB.Flagkk     = Flagkk;
end
SwE.VM         = VM;                %-Filehandle - Mask

SwE.xX         = xX;                %-design structure
SwE.xM         = xM;                %-mask structure

SwE.swd        = pwd;

if isfield(SwE.type,'modified')
  
  SwE.Vis.nCov_vis_g = nCov_vis_g;
  SwE.Vis.nCov_vis = nCov_vis;

  SwE.Gr.nSubj_g   = nSubj_g;
  SwE.Gr.uSubj_g   = uSubj_g;
  
end

%-Save analysis parameters in SwE.mat file
%--------------------------------------------------------------------------
if exist('OCTAVE_VERSION','builtin')
  save('SwE','SwE');
elseif spm_matlab_version_chk('7') >=0
  save('SwE','SwE','-V6');
else
  save('SwE','SwE');
end

%==========================================================================
%- Produce bootstraps and maximum stats 
%==========================================================================

% Produce the random value following the Rademacher distribution
resamplingMatrix = NaN(nScan,WB.nB);
for iS = 1:nSubj
%     resamplingMatrix(iSubj == uSubj(iS),:) = repmat(binornd(1, 0.5, 1, WB.nB), sum(iSubj == uSubj(iS)), 1);
  resamplingMatrix(iSubj == uSubj(iS),:) = repmat(randi([0 1], 1, WB.nB), sum(iSubj == uSubj(iS)), 1);  % BG (08/11/2016): using randi instead of binornd (which is from the stats toolbox)
end
resamplingMatrix(resamplingMatrix == 0) = -1;

% load original score
if isMat
  originalScore = score;
  clear score;
else
  originalScore = spm_get_data(Vscore, XYZ);
  % # blocks
  blksz  = ceil(mmv);                             %-block size
  nbch   = ceil(S/ blksz);          
end
% variables for results
uncP = ones(1, S); % one because of the original score

str   = sprintf('Parameter estimation\nBootstraping');
spm_progress_bar('Init',100,str,'');

fprintf('\n')
for b = 1:WB.nB
  tic
  str   = sprintf('Parameter estimation\nBootstrap # %i', b);
  spm_progress_bar('Set','xlabel', str)

  % activated voxels for cluster-wise inference
  if (SwE.WB.clusterWise == 1)
    activatedVoxels = false(1,S);
    if (SwE.WB.stat == 'T')
      activatedVoxelsNeg = false(1,S);
    end
  end
  if ~isMat
    for bch = 1:nbch                     %-loop over blocks
      blksz  = ceil(mmv);                             %-block size
      if bch ~= nbch
        index = (1+(bch-1)*blksz) : (bch * blksz);
        count = (bch-1)*blksz;
      else
        index = (1+(bch-1)*blksz) : size(XYZ,2);
        blksz = length(index);
      end
      %-Print progress information in command window
      %------------------------------------------------------------------
      str = sprintf('Bootstrap # %i  Block %i/%i', b, bch, nbch);
      
      if  bch == 1
        str2 = '';
      else
        str2 = repmat(sprintf('\b'),1,43);
      end
      fprintf('%s%-40s: %1s',str2,str,' ');
      
      Y_b = spm_get_data(VYWB, XYZ(:,index),false) + ...
        spm_get_data(VResWB, XYZ(:,index),false) .* repmat(resamplingMatrix(:,b),1,blksz);
      
      beta  = pX * Y_b;                     %-Parameter estimates
      if WB.RSwE == 0
        res=swe_fit(SwE, Y_b, xX.X, corr, beta, SwE.SS);
      else 
        res=swe_fit(SwE, Y_b, tmpR2, corr, beta, SwE.SS);
      end
      
      clear Y_b
      
      %-Estimation of the data variance-covariance components (modified SwE)
      %-SwE estimation (classic version)
      %--------------------------------------------------------------
      if isfield(SwE.type,'modified')
        Cov_vis=zeros(nCov_vis,blksz);
        for i = Ind_Cov_vis_diag
          Cov_vis(i,:) = mean(res(Flagk(i,:),:).^2);
        end
        for i = Ind_Cov_vis_off_diag
          if any(Flagk(i,:))
            Cov_vis(i,:)= sum(res(Flagk(i,:),:).*res(Flagkk(i,:),:), 1).*...
              sqrt(Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,1)),:).*...
              Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,2)),:)./...
              sum(res(Flagk(i,:),:).^2, 1)./...
              sum(res(Flagkk(i,:),:).^2, 1));
          end
        end
        %NaN may be produced in cov. estimation when one correspondant
        %variance are = 0, so set them to 0
        Cov_vis(isnan(Cov_vis))=0;
        %need to check if the eigenvalues of Cov_vis matrices are >=0
        for g = 1:SwE.Gr.nGr
          for iVox = 1:blksz
            tmp = zeros(nVis_g(g));
            tmp(tril(ones(nVis_g(g)))==1) = Cov_vis(iGr_Cov_vis_g==g,iVox);
            tmp = tmp + tmp' - diag(diag(tmp));
            [V D] = eig(tmp);
            if any (diag(D)<0) %Bug corrected (BG - 19/09/13)
              D(D<0) = 0;
              tmp = V * D * V';
              Cov_vis(iGr_Cov_vis_g==g,iVox) = tmp(tril(ones(nVis_g(g)))==1); %Bug corrected (BG - 19/09/13)
            end
          end
        end
        cCovBc = weightR * Cov_vis;
      else
        cCovBc = 0;
        for i = 1:nSubj
          Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
            (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
          cCovBc = cCovBc + Cov_beta_i_tmp;
        end
      end
      
      % compute the score
      if (SwE.WB.stat == 'T')
          
        score = (conWB * beta) ./ sqrt(cCovBc);
        clear beta
        
      else

        cBeta = conWB * beta;
        clear beta
        score = zeros(1, blksz);
        for iVox = 1:blksz
          cCovBc_vox = zeros(nSizeCon);
          cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
          cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
          score(iVox) = cBeta(:,iVox)' / cCovBc_vox * cBeta(:,iVox);
        end
        score = score / rankCon;
        
      end
      
      % hypothesis test
      if (WB.clusterWise == 1)
        [~, activatedVoxels(index)]=swe_hyptest(SwE, score, blksz, edf, cCovBc, Cov_vis, dofMat);
        clear cCovBc
      end
      uncP(index) = uncP(index) + (score >= originalScore(index));
          
      maxScore(b+1) = max(maxScore(b+1), max(score));
      if (SwE.WB.stat == 'T')
         minScore(b+1) = min(score);
      end
      
    end % (bch)
  else
    
    %-Print progress information in command window
    %------------------------------------------------------------------
    str = sprintf('Bootstrap # %i', b);
    
    fprintf('%-40s: %1s',str,' ');
    
    Y_b = YWB + resWB .* repmat(resamplingMatrix(:,b),1,S);
    
    beta  = pX * Y_b;                     %-Parameter estimates
    if WB.RSwE == 0
      res=swe_fit(SwE, Y_b, xX.X, corr, beta, SwE.SS);
    else 
      res=swe_fit(SwE, Y_b, tmpR2, corr, beta, SwE.SS);
    end

    clear Y_b
    
    %-Estimation of the data variance-covariance components (modified SwE)
    %-SwE estimation (classic version)
    %--------------------------------------------------------------
    if isfield(SwE.type,'modified')
      Cov_vis=zeros(nCov_vis,S);
      for i = Ind_Cov_vis_diag
        Cov_vis(i,:) = mean(res(Flagk(i,:),:).^2);
      end
      for i = Ind_Cov_vis_off_diag
        if any(Flagk(i,:))
          Cov_vis(i,:)= sum(res(Flagk(i,:),:).*res(Flagkk(i,:),:), 1).*...
            sqrt(Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,1)),:).*...
            Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,2)),:)./...
            sum(res(Flagk(i,:),:).^2, 1)./...
            sum(res(Flagkk(i,:),:).^2, 1));
        end
      end
      %NaN may be produced in cov. estimation when one correspondant
      %variance are = 0, so set them to 0
      Cov_vis(isnan(Cov_vis))=0;
      %need to check if the eigenvalues of Cov_vis matrices are >=0
      for g = 1:SwE.Gr.nGr
        for iVox = 1:S
          tmp = zeros(nVis_g(g));
          tmp(tril(ones(nVis_g(g)))==1) = Cov_vis(iGr_Cov_vis_g==g,iVox);
          tmp = tmp + tmp' - diag(diag(tmp));
          [V D] = eig(tmp);
          if any (diag(D)<0) %Bug corrected (BG - 19/09/13)
            D(D<0) = 0;
            tmp = V * D * V';
            Cov_vis(iGr_Cov_vis_g==g,iVox) = tmp(tril(ones(nVis_g(g)))==1); %Bug corrected (BG - 19/09/13)
          end
        end
      end
      cCovBc = weightR * Cov_vis;
    else
      cCovBc = 0;
      for i = 1:nSubj
        Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
          (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
        cCovBc = cCovBc + Cov_beta_i_tmp;
      end
    end
    
    % compute the score
    if (SwE.WB.stat == 'T')

      score = (conWB * beta) ./ sqrt(cCovBc);
      clear beta
      
    else
      
      cBeta = conWB * beta;
      clear beta
      score = zeros(1, S);
      for iVox = 1:S
        cCovBc_vox = zeros(nSizeCon);
        cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
        cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
        score(iVox) = cBeta(:,iVox)' / cCovBc_vox * cBeta(:,iVox);
      end
      score = score / rankCon;     
    end
    
    % hypothesis test
    if (WB.clusterWise == 1)
        [~, activatedVoxels, activatedVoxelsNeg]=swe_hyptest(SwE, score, S, edf, cCovBc, Cov_vis, dofMat);
        clear cCovBc
    end
    uncP = uncP + (score >= originalScore); 
    
    maxScore(b+1) = max(score);
    if (SwE.WB.stat == 'T')
      minScore(b+1) = min(score);
    end
        
    
  end
  
  % compute the max cluster size if needed (so many ways this can be
  % done... Not sure this solution is the best)
  if (WB.clusterWise == 1)
    if ~isMat || isfield(SwE.WB.clusterInfo, 'Vxyz')      
      LocActivatedVoxels = XYZ(:,activatedVoxels);
      clusterAssignment = spm_clusters(LocActivatedVoxels);
    else %surface data
      LocActivatedVoxels = false(nVox,1);
      LocActivatedVoxels(Cm) = activatedVoxels;
      clusterAssignment = spm_mesh_clusters(faces,LocActivatedVoxels);
      clusterAssignment = clusterAssignment(LocActivatedVoxels)';
    end
    nCluster     = max(clusterAssignment);
    clusterSize = histc(clusterAssignment,1:nCluster);
    if isempty(clusterSize)
      maxClusterSize(b+1) = 0;
    else
      maxClusterSize(b+1) = max(clusterSize);
    end
    if (WB.stat == 'T')
      if ~isMat || isfield(SwE.WB.clusterInfo, 'Vxyz') 
        LocActivatedVoxelsNeg = XYZ(:,activatedVoxelsNeg);
        clusterAssignmentNeg = spm_clusters(LocActivatedVoxelsNeg);
      else %surface data
        LocActivatedVoxelsNeg = false(nVox,1);
        LocActivatedVoxelsNeg(Cm) = activatedVoxelsNeg;
        clusterAssignmentNeg = spm_mesh_clusters(faces,LocActivatedVoxelsNeg);
        clusterAssignmentNeg = clusterAssignmentNeg(LocActivatedVoxelsNeg)';
        if isnan(clusterAssignmentNeg)
          clusterAssignmentNeg = [];
        end
      end
      nClusterNeg     = max(clusterAssignmentNeg);
      clusterSizeNeg = histc(clusterAssignmentNeg,1:nClusterNeg);
      if isempty(clusterSizeNeg)
        maxClusterSizeNeg(b+1) = 0;
      else
        maxClusterSizeNeg(b+1) = max(clusterSizeNeg);
      end
    end
  end
  toc
  spm_progress_bar('Set',100 * b / WB.nB);
end

%-Save analysis original max min in SwE structure
%--------------------------------------------------------------------------
SwE.WB.maxScore = maxScore;
if (WB.clusterWise == 1)
    SwE.WB.clusterInfo.maxClusterSize = maxClusterSize;
end
if (WB.stat == 'T')
    SwE.WB.minScore = minScore;
    if (WB.clusterWise == 1)
        SwE.WB.clusterInfo.maxClusterSizeNeg = maxClusterSizeNeg;
    end
end

%==========================================================================
%- produce results images
%==========================================================================
if isMat
  uncP = uncP / (WB.nB + 1);
  uncP_pos = nan(1, nVox);
  uncP_pos(:,Cm) = uncP;
  VlP_wb_pos = -log10(uncP);
  save(sprintf('swe_vox_%cstat_lp-WB_c%02d%s', WB.stat, 1, file_ext), 'VlP_wb_pos');
  clear VlP_wb_pos
  
  if WB.stat == 'T'
    uncP_neg = 1 + 1/(WB.nB + 1) - uncP_pos;
    VlP_wb_neg = -log10(uncP_neg);
    save(sprintf('swe_vox_%cstat_lp-WB_c%02d%s', WB.stat, 2, file_ext), 'VlP_wb_neg');
    clear VlP_wb_neg
    
    VZ_wb = swe_invNcdf(1 - uncP);
    save(sprintf('swe_vox_z%cstat-WB_c%02d%s', WB.stat, 1, file_ext), 'VZ_wb');
    clear VZ_wb
    
  else
      
    VX_wb = spm_invXcdf(1 - uncP,1);
    save(sprintf('swe_vox_x%cstat-WB_c%02d%s', WB.stat, 1, file_ext), 'VX_wb');
    clear VX_wb
    
  end
  
  %
  % - write out lP_FWE+ and lP_FWE- ;
  %
  tol = 1e-4;	% Tolerance for comparing real numbers
  
  FWERP = ones(1, S); % 1 because the original maxScore is always > original Score
  
  for b = 1:WB.nB
    %-FWER-corrected p is proportion of randomisation greater or
    % equal to statistic.
    %-Use a > b -tol rather than a >= b to avoid comparing
    % two reals for equality.
    FWERP = FWERP + (maxScore(b+1) > originalScore - tol);
  end
  FWERP = FWERP / (WB.nB + 1);
  fwerP_pos = nan(1, nVox);
  fwerP_pos(:,Cm) = FWERP;
  VlP_wb_FWE_pos = -log10(fwerP_pos);
  save(sprintf('swe_vox_%cstat_lpFWE-WB_c%02d%s',WB.stat,1,file_ext), 'VlP_wb_FWE_pos');
  clear VlP_wb_FWE_pos fwerP_pos FWERP
  
  
  if WB.stat == 'T'
    FWERPNeg = ones(1, S); % 1 because the original maxScore is always > original Score
    
    for b = 1:WB.nB
      %-FWER-corrected p is proportion of randomisation greater or
      % equal to statistic.
      %-Use a > b -tol rather than a >= b to avoid comparing
      % two reals for equality.
      FWERPNeg = FWERPNeg + (minScore(b+1) < originalScore + tol);
    end
    FWERPNeg = FWERPNeg / (WB.nB + 1);
    fwerP_neg = nan(1, nVox);
    fwerP_neg(:,Cm) = FWERPNeg;
    VlP_wb_FWE_neg = -log10(fwerP_neg);
    save(sprintf('swe_vox_%cstat_lpFWE-WB_c%02d%s',WB.stat,2,file_ext), 'VlP_wb_FWE_neg');
    clear VlP_wb_FWE_neg fwerP_neg
  end
  
  %
  % - write out lP_FDR+ and lP_FDR- images;
  %
  try
    fdrP = spm_P_FDR(uncP);
  catch
    fdrP = spm_P_FDR(uncP,[],'P',[],sort(uncP)');
  end
  fdrP_pos = nan(1, nVox);
  fdrP_pos(:,Cm) = fdrP;
  VlP_wb_FDR_pos = -log10(fdrP_pos);
  save(sprintf('swe_vox_%cstat_lpFDR-WB_c%02d%s',WB.stat,1,file_ext), 'VlP_wb_FDR_pos');
  clear VlP_wb_FDR_pos fdrP_pos fdrP
  
  if WB.stat =='T'
    try
      fdrP = spm_P_FDR(1 + 1/(WB.nB + 1) - uncP);
    catch
      fdrP = spm_P_FDR(1 + 1/(WB.nB + 1) - uncP,[],'P',[],sort(1 + 1/(WB.nB + 1) - uncP)');
    end
    fdrP_neg = nan(1, nVox);
    fdrP_neg(:,Cm) = fdrP;
    VlP_wb_FDR_neg = -log10(fdrP_neg);
    save(sprintf('swe_vox_%cstat_lpFDR-WB_c%02d%s',WB.stat,2,file_ext), 'VlP_wb_FDR_neg');
    clear VlP_wb_FDR_neg fdrP_neg fdrP
  end
  
  if WB.clusterWise == 1
    % Not sure what to output. So might be changed later.
    % For now, -log(p_{cluster-wise FWER}) image with nan for non-surviving
    % voxels after the thresholding of the original data
    
    clusterFwerP_pos_perCluster = ones(1, SwE.WB.clusterInfo.nCluster); % 1 because the original maxScore is always > original Score
    if (~isempty(SwE.WB.clusterInfo.clusterSize))
      for b = 1:WB.nB
        clusterFwerP_pos_perCluster = clusterFwerP_pos_perCluster + (maxClusterSize(b+1) >= SwE.WB.clusterInfo.clusterSize);
      end
      clusterFwerP_pos_perCluster = clusterFwerP_pos_perCluster / (WB.nB + 1);
    end
    
    clusterFwerP_pos_perElement = nan(1, nVox);
    if ~isMat || isfield(SwE.WB.clusterInfo, 'Vxyz')      
      tmp = find(Cm);
      tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxels,2));
      for iC = 1:SwE.WB.clusterInfo.nCluster
        tmp3(SwE.WB.clusterInfo.clusterAssignment == iC) = clusterFwerP_pos_perCluster(iC);
      end
      clusterFwerP_pos_perElement(tmp(activatedVoxels)) = tmp3;
    else
      tmp3 = zeros(1, sum(SwE.WB.clusterInfo.LocActivatedVoxels));
      for iC = 1:SwE.WB.clusterInfo.nCluster
        tmp3(SwE.WB.clusterInfo.clusterAssignment == iC) = clusterFwerP_pos_perCluster(iC);
      end
      clusterFwerP_pos_perElement(SwE.WB.clusterInfo.LocActivatedVoxels) = tmp3;
    end
    VlP_wb_clusterFWE_pos  = -log10(clusterFwerP_pos_perElement);
    save(sprintf('swe_clustere_%cstat_lpFWE-WB_c%02d%s',WB.stat,1,file_ext), 'VlP_wb_clusterFWE_pos');
    
    if WB.stat =='T'
      
      clusterFwerP_neg_perCluster = ones(1, SwE.WB.clusterInfo.nClusterNeg); % 1 because the original maxScore is always > original Score
      if (~isempty(SwE.WB.clusterInfo.clusterSizeNeg))
        for b = 1:WB.nB
          clusterFwerP_neg_perCluster = clusterFwerP_neg_perCluster + (maxClusterSizeNeg(b+1) >= SwE.WB.clusterInfo.clusterSizeNeg);
        end
        clusterFwerP_neg_perCluster = clusterFwerP_neg_perCluster / (WB.nB + 1);
      end
      
      clusterFwerP_neg_perElement = nan(1, nVox);
      if ~isMat || isfield(SwE.WB.clusterInfo, 'Vxyz')
        tmp = find(Cm);
        tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxelsNeg,2));
        for iC = 1:SwE.WB.clusterInfo.nClusterNeg
          tmp3(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = clusterFwerP_neg_perCluster(iC);
        end
      clusterFwerP_neg_perElement(tmp(activatedVoxelsNeg)) = tmp3;
      else
        tmp3 = zeros(1, sum(SwE.WB.clusterInfo.LocActivatedVoxelsNeg));
        for iC = 1:SwE.WB.clusterInfo.nClusterNeg
          tmp3(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = clusterFwerP_neg_perCluster(iC);
        end
        clusterFwerP_neg_perElement(SwE.WB.clusterInfo.LocActivatedVoxelsNeg) = tmp3;
      end
      VlP_wb_clusterFWE_neg  = -log10(clusterFwerP_neg_perElement);
      save(sprintf('swe_clustere_%cstat_lpFWE-WB_c%02d%s',WB.stat,2,file_ext), 'VlP_wb_clusterFWE_neg');
    end
  end
else
  
  Q = cumprod([1,SwE.xVol.DIM(1:2)'])*XYZ - ...
    sum(cumprod(SwE.xVol.DIM(1:2)'));
  %
  % - write out lP+ and lP- images;
  %
  uncP = uncP / (WB.nB + 1);
  tmp= nan(SwE.xVol.DIM');
  tmp(Q) = -log10(uncP);
  spm_write_vol(VlP_wb_pos, tmp);
  
  % If it's F, write out an X map.
  stat = nan(SwE.xVol.DIM');
  if WB.stat == 'F'
      stat(Q) = spm_invXcdf(1 - uncP,1);
      spm_write_vol(VcScore_wb_pos, stat);
  end
  
  % If it's T, write out a Z map.
  if WB.stat == 'T'
      
    % Positive map.
    stat(Q) = swe_invNcdf(1 - uncP);
    spm_write_vol(VcScore_wb_pos, stat);

    % T is two tailed so we need a negative map as well.
    tmp(Q) = -log10(1 + 1/(WB.nB + 1) - uncP);
    spm_write_vol(VlP_wb_neg, tmp);
  end
  
  %
  % - write out lP_FWE+ and lP_FWE- images;
  %
  tol = 1e-4;	% Tolerance for comparing real numbers
  
  FWERP = ones(1, S); % 1 because the original maxScore is always > original Score
  
  for b = 1:WB.nB
    %-FWER-corrected p is proportion of randomisation greater or
    % equal to statistic.
    %-Use a > b -tol rather than a >= b to avoid comparing
    % two reals for equality.
    FWERP = FWERP + (maxScore(b+1) > originalScore - tol);
  end
  FWERP = FWERP / (WB.nB + 1);
  tmp(Q) = -log10(FWERP);
  spm_write_vol(VlP_wb_FWE_pos, tmp);
  
  if WB.stat == 'T'
    FWERPNeg = ones(1, S); % 1 because the original maxScore is always > original Score
    
    for b = 1:WB.nB
      %-FWER-corrected p is proportion of randomisation greater or
      % equal to statistic.
      %-Use a > b -tol rather than a >= b to avoid comparing
      % two reals for equality.
      FWERPNeg = FWERPNeg + (minScore(b+1) < originalScore + tol);
    end
    FWERPNeg = FWERPNeg / (WB.nB + 1);
    tmp(Q) = -log10(FWERPNeg);
    spm_write_vol(VlP_wb_FWE_neg, tmp);
  end
  
  %
  % - write out lP_FDR+ and lP_FDR- images;
  %
  try
    tmp(Q) = -log10(spm_P_FDR(uncP));
  catch
    tmp(Q) = -log10(spm_P_FDR(uncP,[],'P',[],sort(uncP)'));
  end
  spm_write_vol(VlP_wb_FDR_pos, tmp);
  
  if WB.stat =='T'
    try
      tmp(Q) = -log10(spm_P_FDR(1 + 1/(WB.nB + 1) - uncP));
    catch
      tmp(Q) = -log10(spm_P_FDR(1 + 1/(WB.nB + 1) - uncP,[],'P',[],sort(1 + 1/(WB.nB + 1) - uncP)'));
    end
    spm_write_vol(VlP_wb_FDR_neg, tmp);
  end
  
  if WB.clusterWise == 1
    % Not sure what to output. So might be changed later.
    % For now, -log(p_{cluster-wise FWER}) image with nan for non-surviving
    % voxels after the thresholding of the original data
    Q = cumprod([1,SwE.xVol.DIM(1:2)']) * SwE.WB.clusterInfo.LocActivatedVoxels - ...
      sum(cumprod(SwE.xVol.DIM(1:2)'));
    tmp= nan(SwE.xVol.DIM');
    
    clusterFwerP_pos_perCluster = ones(1, SwE.WB.clusterInfo.nCluster); % 1 because the original maxScore is always > original Score
    if (~isempty(SwE.WB.clusterInfo.clusterSize))
      for b = 1:WB.nB
        clusterFwerP_pos_perCluster = clusterFwerP_pos_perCluster + (maxClusterSize(b+1) >= SwE.WB.clusterInfo.clusterSize);
      end
      clusterFwerP_pos_perCluster = clusterFwerP_pos_perCluster / (WB.nB + 1);
    end
    tmp2 = -log10(clusterFwerP_pos_perCluster);
    
    tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxels,2));
    for iC = 1:SwE.WB.clusterInfo.nCluster
      tmp3(SwE.WB.clusterInfo.clusterAssignment == iC) = tmp2(iC);
    end
    tmp(Q) = tmp3;
    spm_write_vol(VlP_wb_clusterFWE_pos, tmp);
    if WB.stat =='T'
      Q = cumprod([1,SwE.xVol.DIM(1:2)']) * SwE.WB.clusterInfo.LocActivatedVoxelsNeg - ...
        sum(cumprod(SwE.xVol.DIM(1:2)'));
      tmp= nan(SwE.xVol.DIM');
      
      clusterFwerP_neg_perCluster = ones(1, SwE.WB.clusterInfo.nClusterNeg); % 1 because the original maxScore is always > original Score
      if (~isempty(SwE.WB.clusterInfo.clusterSizeNeg))
        for b = 1:WB.nB
          clusterFwerP_neg_perCluster = clusterFwerP_neg_perCluster + (maxClusterSizeNeg(b+1) >= SwE.WB.clusterInfo.clusterSizeNeg);
        end
        clusterFwerP_neg_perCluster = clusterFwerP_neg_perCluster / (WB.nB + 1);
      end
      tmp2 = -log10(clusterFwerP_neg_perCluster);
      
      tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxelsNeg, 2));
      for iC = 1:SwE.WB.clusterInfo.nClusterNeg
        tmp3(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = tmp2(iC);
      end
      tmp(Q) = tmp3;
      spm_write_vol(VlP_wb_clusterFWE_neg, tmp);
    end
  end
end
      
%==========================================================================
%- E N D: Cleanup GUI
%==========================================================================

if ~isMat
    % Remove residual and Y images now we are done with them:
    files = {'^swe_vox_resid_y.{4}\..{3}$','^swe_vox_fit_y.{4}\..{3}$'};
    for i = 1:numel(files)
      j = cellstr(spm_select('FPList',SwE.swd,files{i}));
      for k = 1:numel(j)
        spm_unlink(j{k});
      end
    end
end

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#
%spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
fprintf('...use the saved images for assessment\n\n')

end

%-Mask out voxels where data is constant in at least one separable
% matrix design either in a visit category or within-subject (BG - 27/05/2016)
function [Cm,Y,CrS]=swe_mask_seperable(SwE, Cm, Y, iGr_dof)
    
      % Setup
      nGr_dof = length(unique(iGr_dof));
      nGr = SwE.Gr.nGr;
      iGr = SwE.Gr.iGr;
      uGr = SwE.Gr.uGr;
      iVis = SwE.Vis.iVis;
      iSubj = SwE.Subj.iSubj;
      nVis_g = SwE.Vis.nVis_g;
      uVis_g = SwE.Vis.uVis_g;
      
      for g = 1:nGr_dof % first look data for each separable matrix design
        if sum(iGr_dof'==g) > 1 % do not look for cases where the separable matrix design is only one row (BG - 05/08/2016)
          Cm(Cm) = any(abs(diff(Y(iGr_dof'==g,Cm),1)) > eps, 1); % mask constant data within separable matrix design g (added by BG on 29/08/16)
          if isfield(SwE.type,'modified') % added by BG on 29/08/16
            for g2 = 1:nGr % then look data for each "homogeneous" group
              % check if the data is contant over subject for each visit category
              for k = 1:nVis_g(g2)
                if sum(iGr_dof'==g & iGr == uGr(g2) & iVis == uVis_g{g2}(k)) > 1 % do not look for cases when the data is only one row (BG - 05/08/2016)
                  Cm(Cm) = any(abs(diff(Y(iGr_dof'==g & iGr == uGr(g2) & iVis == uVis_g{g2}(k) ,Cm),1)) > eps, 1);
                  for kk = k:nVis_g(g2)
                    if k ~= kk
                      % extract the list of subject with both visit k and kk
                      subjList = intersect(iSubj(iGr_dof'==g & iGr == uGr(g2) & iVis == uVis_g{g2}(k)), iSubj(iGr_dof'==g & iGr == uGr(g2) & iVis == uVis_g{g2}(kk)));
                      % look if some difference are observed within subject
                      if ~isempty(subjList)
                        diffVis = Cm(Cm) == 0;
                        for i = 1:length(subjList)
                          diffVis = diffVis | (abs(Y(iSubj == subjList(i) & iVis == uVis_g{g2}(k), Cm) - Y(iSubj == subjList(i) & iVis == uVis_g{g2}(kk), Cm)) > eps);
                        end
                        Cm(Cm) = diffVis;
                      end
                    end
                  end
                end
              end
            end
          end
        end
      end
      
      Y      = Y(:,Cm);                          %-Data within mask
      CrS    = sum(Cm);                          %-# current voxels

end

% This function performs a hypothesis test using the threshold given as the
% primary threshold in the SwE cluster info. If this is not available it
% returns unthresholded p values only.
function [p, activatedVoxels, activatedVoxelsNeg]=swe_hyptest(SwE, score, matSize, edf, cCovBc, Cov_vis, dofMat, varargin)

      % setup
      p = zeros(1, matSize);
      nGr = length(unique(SwE.Gr.iGr));
      nSizeCon = size(SwE.WB.con,1);
      rankCon = rank(SwE.WB.con);

      if nSizeCon == 1
          Wg_2 = SwE.WB.Wg{1};
	      Wg_3 = SwE.WB.Wg{1};
      else
	      Wg_2 = SwE.WB.Wg{2};
	      Wg_3 = SwE.WB.Wg{3};
      end

      if isfield(SwE.type,'modified')
         dof_type = SwE.type.modified.dof_mo;
      else
         dof_type = SwE.type.classic.dof_cl;        
      end

      % Convert P values.
	  switch dof_type
	     case 1
	        error('degrees of freedom type still not implemented for the WB')

	     case 2
	        CovcCovBc = 0;
	        for g = 1:nGr
	          CovcCovBc = CovcCovBc + Wg_2{g} * swe_vechCovVechV(Cov_vis(SwE.WB.iGr_Cov_vis_g==g,:), dofMat{g}, 1);
	        end
	        if (SwE.WB.stat == 'T')
	           edf = 2 * cCovBc.^2 ./ CovcCovBc - 2;
	        else
	           edf = 2 * (sum(swe_duplication_matrix(nSizeCon), 1) * cCovBc).^2 ./ CovcCovBc - 2;
	        end

	     case 3
	        CovcCovBc = 0;
	        for g = 1:nGr
	          CovcCovBc = CovcCovBc + Wg_3{g} * swe_vechCovVechV(Cov_vis(SwE.WB.iGr_Cov_vis_g==g,:), dofMat{g}, 2);
	        end
	        if (SwE.WB.stat == 'T')
    	       edf = 2 * cCovBc.^2 ./ CovcCovBc;
    	    else
			   tmp = eye(nSizeCon);
               edf = (sum(swe_duplication_matrix(nSizeCon), 1) * cCovBc.^2 +...
                     (tmp(:)' * swe_duplication_matrix(nSizeCon) * cCovBc).^2) ./ CovcCovBc;
    	    end
	   end

	   % P values and activated voxels (if clusterwise).
	   if (SwE.WB.stat == 'T')
    	  p  = spm_Tcdf(score, edf);
          
          if SwE.WB.clusterWise~=0
              if nargin <=7
                % We may wish to just record the activated voxels. 
                activatedVoxels = p > (1-SwE.WB.clusterInfo.primaryThreshold);
                activatedVoxelsNeg = p < (SwE.WB.clusterInfo.primaryThreshold);
              else
                % Or we may wish to add the activatedVoxels to a pre-existing list.
                activatedVoxels = [varargin{1}, p > (1-SwE.WB.clusterInfo.primaryThreshold)];
                activatedVoxelsNeg = [varargin{2}, p < (SwE.WB.clusterInfo.primaryThreshold)];
              end
          end

	   else
	   	  scoreTmp = (edf-rankCon+1) ./ edf .* score;
	      scoreTmp(scoreTmp < 0 ) = 0;
	      if dof_type == 0
	        p(scoreTmp>0) = 1-betainc((edf-rankCon+1)./(edf-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf-rankCon+1)/2, rankCon/2);
	      else
	        p(scoreTmp>0) = 1-betainc((edf(scoreTmp>0)-rankCon+1)./(edf(scoreTmp>0)-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf(scoreTmp>0)-rankCon+1)/2, rankCon/2);
	        p(scoreTmp == 0) = 0;
          end

          if SwE.WB.clusterWise~=0
              if nargin<=7
                  activatedVoxels = p > (1-SwE.WB.clusterInfo.primaryThreshold);
              else
                  activatedVoxels = [varargin{1}, p > (1-SwE.WB.clusterInfo.primaryThreshold)];
              end
          end
          
          activatedVoxelsNeg = NaN;
          
	   end

end

% This function performs the users requested residual corrections and
% calculates tmpR2 (the adjusted xX.X).
function [corr, tmpR2] = swe_resid_corr(SwE, restric, ss, pX, varargin)

    xX = SwE.xX;
    [nScan, nBeta] = size(xX.X);
    conWB = SwE.WB.con;
    iSubj = SwE.Subj.iSubj;
    nSubj = SwE.Subj.nSubj;
    uSubj = SwE.Subj.uSubj;
    rankCon = rank(SwE.WB.con);

    % This is to prevent tmpR2 being overwritten.
    if nargin <= 4
        tmpR2 = false;
    else
        tmpR2 = varargin{1};
    end
    
    if restric == 1
      tmpR = (xX.X' * xX.X) \ conWB';
      tmpR = tmpR / (conWB * tmpR);
      tmpR2 = xX.X * (eye(nBeta) - tmpR * conWB);
      Hat = xX.X * (pX - tmpR * conWB * pX); % Restricted Hat matrix
    else
      Hat = xX.X*(pX); % Hat matrix
    end

    switch ss
        case 0
          corr = ones(nScan,1);
        case 1
          if WB.RSwE == 1
            corr  = repmat(sqrt(nScan/(nScan - nBeta + rankCon)),nScan,1); % residual correction (type 1)
          else
            corr  = repmat(sqrt(nScan/(nScan-nBeta)),nScan,1); 
          end
        case 2
          corr  = (1-diag(Hat)).^(-0.5); % residual correction (type 2)
        case 3
          corr  = (1-diag(Hat)).^(-1); % residual correction (type 3)
        case 4
          corr =cell(nSubj,1);
          I_Hat = eye(nScan) - Hat;
          for i = 1:nSubj
            tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
            tmp = (tmp + tmp')/2;
            [tmpV, tmpE] = eig(tmp);
            corr{i} = tmpV * diag(1./sqrt(diag(tmpE))) * tmpV';
          end
        case 5
          corr  = cell(nSubj,1);
          I_Hat = eye(nScan) - Hat;
          for i = 1:nSubj
            tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
            tmp = (tmp + tmp')/2;
            corr{i} = inv(tmp);
          end
    end
end

% This function obtains Y estimates and residuals from fitting data.
function [res, Y_est]=swe_fit(SwE, Y, crctX, corr, beta, ss)
    
    Y_est = crctX * beta;
    if ss >= 4 % SC2 or SC3
      res = zeros(size(Y));
      for i = 1:SwE.Subj.nSubj
        res(SwE.Subj.iSubj==SwE.Subj.uSubj(i),:) = corr{i} *...
            (Y(SwE.Subj.iSubj==SwE.Subj.uSubj(i),:)-...
            Y_est(SwE.Subj.iSubj==SwE.Subj.uSubj(i),:));
      end
    else
      res  = diag(corr) * (Y-Y_est);
    end
    
end