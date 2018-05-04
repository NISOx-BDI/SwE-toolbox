function swe_cp_WB(SwE)

%-Say hello
%--------------------------------------------------------------------------
Finter = spm('CreateIntWin','off');
set(Finter,'name','SwE estimation');
set(Finter,'vis','on')

%-Change to SwE.swd if specified
%--------------------------------------------------------------------------
try
  cd(SwE.swd);
catch %#ok<*CTCH>
  SwE.swd = pwd;
end

%-Ensure data are assigned
%--------------------------------------------------------------------------
try
  SwE.xY.VY;
catch
  spm('alert!','Please assign data to this design', mfilename);
  spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
  return
end

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
if ~isempty(spm_select('List',SwE.swd,'^mask\..{3}$'))
  
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
  
% small sample correction 
if WB.RWB == 1
  tmpR = (xX.X' * xX.X) \ conWB';
  tmpR = tmpR / (conWB * tmpR);
  tmpR2 = xX.X * (eye(nBeta) - tmpR * conWB);
  HatR = xX.X * (pX - tmpR * conWB * pX);
  
  switch WB.SS
    case 0 
      corrWB = ones(nScan,1);
    case 1
      corrWB  = repmat(sqrt(nScan/(nScan - nBeta + rankCon)),nScan,1); % residual correction (type 1) 
    case 2
      corrWB  = (1-diag(HatR)).^(-0.5); % residual correction (type 2)
    case 3
      corrWB  = (1-diag(HatR)).^(-1); % residual correction (type 3)
    case 4
      corrWB =cell(nSubj,1);
      I_Hat = eye(nScan) - HatR;
      for i = 1:nSubj
        tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
        tmp = (tmp + tmp')/2;
        [tmpV, tmpE] = eig(tmp);
        corrWB{i} = tmpV * diag(1./sqrt(diag(tmpE))) * tmpV';
      end
      clear I_Hat tmp
    case 5
      corrWB  = cell(nSubj,1);
      I_Hat = eye(nScan) - HatR;
      for i = 1:nSubj
          tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
          tmp = (tmp + tmp')/2;
          corrWB{i} = inv(tmp); 
      end
      clear I_Hat tmp
  end
else 
  Hat = xX.X*(pX); % Hat matrix
  switch WB.SS
    case 0
      corrWB = ones(nScan,1);
    case 1
      corrWB  = repmat(sqrt(nScan/(nScan-nBeta)),nScan,1); % residual correction (type 1)
    case 2
      corrWB  = (1-diag(Hat)).^(-0.5); % residual correction (type 2)
    case 3
      corrWB  = (1-diag(Hat)).^(-1); % residual correction (type 3)
    case 4
      corrWB  = cell(nSubj,1);
      I_Hat = eye(nScan) - Hat;
      for i = 1:nSubj
        tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
        tmp = (tmp + tmp')/2;
        [tmpV, tmpE] = eig(tmp);
        corrWB{i} = tmpV * diag(1./sqrt(diag(tmpE))) * tmpV';
      end
      clear I_Hat tmp
    case 5
      corrWB  = cell(nSubj,1);
      I_Hat = eye(nScan) - Hat;
      for i = 1:nSubj
        tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
        tmp = (tmp + tmp')/2;
        corrWB{i} = inv(tmp);
      end
      clear I_Hat tmp
  end
end

if WB.RSwE == 1
  tmpR = (xX.X' * xX.X) \ conWB';
  tmpR = tmpR / (conWB * tmpR);
  tmpR2 = xX.X * (eye(nBeta) - tmpR * conWB);
  HatR = xX.X * (pX - tmpR * conWB * pX);
  
  switch SwE.SS
    case 0
      corr = ones(nScan,1);
    case 1
      corr  = repmat(sqrt(nScan/(nScan - nBeta + rankCon)),nScan,1); % residual correction (type 1)
    case 2
      corr  = (1-diag(HatR)).^(-0.5); % residual correction (type 2)
    case 3
      corr  = (1-diag(HatR)).^(-1); % residual correction (type 3)
    case 4
      corr =cell(nSubj,1);
      I_Hat = eye(nScan) - HatR;
      for i = 1:nSubj
        tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
        tmp = (tmp + tmp')/2;
        [tmpV, tmpE] = eig(tmp);
        corr{i} = tmpV * diag(1./sqrt(diag(tmpE))) * tmpV';
      end
      clear I_Hat tmp
    case 5
      corr  = cell(nSubj,1);
      I_Hat = eye(nScan) - HatR;
      for i = 1:nSubj
        tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
        tmp = (tmp + tmp')/2;
        corr{i} = inv(tmp);
      end
      clear I_Hat tmp
  end
else
  Hat = xX.X*(pX); % Hat matrix
  switch SwE.SS
    case 0
      corr = ones(nScan,1);
    case 1
      corr  = repmat(sqrt(nScan/(nScan-nBeta)),nScan,1); % residual correction (type 1)
    case 2
      corr  = (1-diag(Hat)).^(-0.5); % residual correction (type 2)
    case 3
      corr  = (1-diag(Hat)).^(-1); % residual correction (type 3)
    case 4
      corr  = cell(nSubj,1);
      I_Hat = eye(nScan) - Hat;
      for i = 1:nSubj
        tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
        tmp = (tmp + tmp')/2;
        [tmpV, tmpE] = eig(tmp);
        corr{i} = tmpV * diag(1./sqrt(diag(tmpE))) * tmpV';
      end
      clear I_Hat tmp
    case 5
      corr  = cell(nSubj,1);
      I_Hat = eye(nScan) - Hat;
      for i = 1:nSubj
        tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
        tmp = (tmp + tmp')/2;
        corr{i} = inv(tmp);
      end
      clear I_Hat tmp
  end
end


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
  
end



%-preprocessing for the modified SwE
if isfield(SwE.type,'modified')
  iVis      = SwE.Vis.iVis;
  iGr       = SwE.Gr.iGr;
  uGr       = unique(iGr);
  nGr       = length(uGr);
  
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
  % Weigth giving only the contrasted SwE (WB)
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
   
%-compute the effective dof from each homogeneous group if dof_type
    switch dof_type
      case 1
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
  VM    = struct('fname',  ['mask' file_ext],...
    'dim',    DIM',...
    'dt',     [spm_type('uint8') spm_platform('bigend')],...
    'mat',    M,...
    'pinfo',  [1 0 0]',...
    'descrip','swe_cp_WB:resultant analysis mask');
  VM    = spm_create_vol(VM);
  
  %-Initialise original score image
  %----------------------------------------------------------------------
  
  Vscore = deal(struct(...
    'fname',    [],...
    'dim',      DIM',...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  ''));
  
  Vscore.fname   = sprintf(['score' file_ext]);
  Vscore.descrip = sprintf('score original data');
  Vscore = spm_create_vol(Vscore);
  
  %-Initialise residual images for the resampling
  %----------------------------------------------------------------------
  VResWB(1:nScan) = deal(struct(...
    'fname',    [],...
    'dim',      DIM',...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  ''));
  
  for i = 1:nScan
    VResWB(i).fname   = sprintf(['ResWB_%04d' file_ext], i);
    if WB.RWB == 1
      VResWB(i).descrip = sprintf('adjusted restricted residuals (%04d)', i);
    else
      VResWB(i).descrip = sprintf('adjusted unrestricted residuals (%04d)', i);
    end
  end
  VResWB = spm_create_vol(VResWB);
  fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');    %-#
  
  %-Initialise fitted data images for the resampling
  %----------------------------------------------------------------------
  VYWB(1:nScan) = deal(struct(...
    'fname',    [],...
    'dim',      DIM',...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  ''));
  
  for i = 1:nScan
    VYWB(i).fname   = sprintf(['YfittedWB_%04d' file_ext], i);
    if WB.RWB == 1
      VYWB(i).descrip = sprintf('restricted fitted data  (%04d)', i);
    else
      VResWB(i).descrip = sprintf('unrestricted fitted data (%04d)', i);
    end
  end
  VYWB = spm_create_vol(VYWB);
  
  %-Initialise result images
  %----------------------------------------------------------------------
  VlP_pos = deal(struct(...
    'fname',    ['lP+' file_ext],...
    'dim',      DIM',...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  '-log10(uncor. non-para. P, +ve)'));
  VlP_pos = spm_create_vol(VlP_pos);
  
  VlP_FWE_pos = deal(struct(...
    'fname',    ['lP_FWE+' file_ext],...
    'dim',      DIM',...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  '-log10(FWE-corr. P, +ve)'));
  VlP_FWE_pos = spm_create_vol(VlP_FWE_pos);
  
  VlP_FDR_pos = deal(struct(...
    'fname',    ['lP_FDR+' file_ext],...
    'dim',      DIM',...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  '-log10(FDR-corr. P, +ve)'));
  VlP_FDR_pos=spm_create_vol(VlP_FDR_pos);
  
  if WB.stat=='T'
    VlP_neg = deal(struct(...
      'fname',    ['lP-' file_ext],...
      'dim',      DIM',...
      'dt',       [spm_type('float32') spm_platform('bigend')],...
      'mat',      M,...
      'pinfo',    [1 0 0]',...
      'descrip',  '-log10(uncor. non-para. P, -ve)'));
    VlP_neg=spm_create_vol(VlP_neg);
    
    VlP_FWE_neg = deal(struct(...
      'fname',    ['lP_FWE-' file_ext],...
      'dim',      DIM',...
      'dt',       [spm_type('float32') spm_platform('bigend')],...
      'mat',      M,...
      'pinfo',    [1 0 0]',...
      'descrip',  '-log10(FWE-corr. P, -ve)'));
    VlP_FWE_neg = spm_create_vol(VlP_FWE_neg);
    
    VlP_FDR_neg = deal(struct(...
      'fname',    ['lP_FDR-' file_ext],...
      'dim',      DIM',...
      'dt',       [spm_type('float32') spm_platform('bigend')],...
      'mat',      M,...
      'pinfo',    [1 0 0]',...
      'descrip',  '-log10(FDR-corr. P, -ve)'));
    VlP_FDR_neg=spm_create_vol(VlP_FDR_neg);
  end
  
  if WB.clusterWise == 1
    VlP_clusterFWE_pos = deal(struct(...
      'fname',    ['lP_clusterFWE+' file_ext],...
      'dim',      DIM',...
      'dt',       [spm_type('float32') spm_platform('bigend')],...
      'mat',      M,...
      'pinfo',    [1 0 0]',...
      'descrip',  '-log10(clusterFWE-corr. P, +ve)'));
    VlP_clusterFWE_pos = spm_create_vol(VlP_clusterFWE_pos);
    
    if WB.stat=='T'
      VlP_clusterFWE_neg = deal(struct(...
        'fname',    ['lP_clusterFWE-' file_ext],...
        'dim',      DIM',...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      M,...
        'pinfo',    [1 0 0]',...
        'descrip',  '-log10(clusterFWE-corr. P, -ve)'));
      VlP_clusterFWE_neg = spm_create_vol(VlP_clusterFWE_neg);
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
  
  for z = 1:nbz:zdim                       %-loop over planes (2D or 3D data)
    
    % current plane-specific parameters
    %----------------------------------------------------------------------
    CrPl         = z:min(z+nbz-1,zdim);       %-plane list
    zords        = CrPl(:)*ones(1,xdim*ydim); %-plane Z coordinates
    CrScore      = [];                        %-scores
    CrYWB        = [];                        %-fitted data under H0
    CrResWB      = [];                        %-residuals
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
        Cm(Cm) = spm_get_data(xM.VM(i),j(:,Cm),false) > 0;
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
      clear diffVis
      Y      = Y(:,Cm);                          %-Data within mask
      CrS    = sum(Cm);                          %-# current voxels
      
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
          YWB = tmpR2 * beta;
          if SwE.WB.SS >= 4 % SC2 or SC3
            resWB = zeros(size(Y));
            for i = 1:nSubj
              resWB(iSubj==uSubj(i),:) = corrWB{i} *...
                (Y(iSubj==uSubj(i),:)-YWB(iSubj==uSubj(i),:));
            end
          else
            resWB  = diag(corrWB) * (Y-YWB);
          end
        else
          YWB = xX.X * beta;
          if SwE.WB.SS >= 4 % SC2 or SC3
            resWB = zeros(size(Y));
            for i = 1:nSubj
              resWB(iSubj==uSubj(i),:) = corrWB{i} *...
                (Y(iSubj==uSubj(i),:)-YWB(iSubj==uSubj(i),:));
            end
          else
            resWB  = diag(corrWB) * (Y-YWB);
          end
        end
        
        if WB.RSwE == 1
          if SwE.SS >= 4  % Cluster-wise adjustments
            res = zeros(size(Y));
            for i = 1:nSubj
              res(iSubj==uSubj(i),:) = corr{i} *...
                (Y(iSubj==uSubj(i),:)-tmpR2(iSubj==uSubj(i),:)*beta);
            end
          else
            res = diag(corr) * (Y - tmpR2 * beta); %-Corrected residuals
          end
        else
          if SwE.SS >= 4  % Cluster-wise adjustments
            res = zeros(size(Y));
            for i = 1:nSubj
              res(iSubj==uSubj(i),:) = corr{i} *...
                (Y(iSubj==uSubj(i),:)-xX.X(iSubj==uSubj(i),:)*beta);
            end
          else
            res = diag(corr) * (Y-xX.X*beta); %-Corrected residuals
          end
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
        else
          cCovBc = 0;
          for i = 1:nSubj
            Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
              (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
            cCovBc = cCovBc + Cov_beta_i_tmp;
          end
        end
        
        % compute the score
        if (nSizeCon == 1)
          if isfield(SwE.type,'modified')
            cCovBc = weightR * Cov_vis;
          end
          
          score = (conWB * beta) ./ sqrt(cCovBc);
          
          if (SwE.WB.clusterWise == 1)
            % need to convert score into parametric p-values
            p = zeros(1, CrS);
            switch dof_type
              case 0
                if WB.stat == 'T'
                  if any(score > 0)
                    p(score > 0)  = spm_Tcdf(-score(score>0), edf);
                  end
                  if any(score <= 0)
                    p(score <= 0) = spm_Tcdf(score(score<=0), edf);
                  end
                else
                  p = 2 * spm_Tcdf(-abs(score), edf);
                end
              case 1
                error('degrees of freedom type still not implemented for the WB')
                
              case 2
                CovcCovBc = 0;
                for g = 1:nGr
                  %CovcCovBc = CovcCovBc + Wg{g} * swe_vechCovVechV(Cov_vis(SwE.Vis.iGr_Cov_vis_g==g,:), dofMat{g}, 1);
                  CovcCovBc = CovcCovBc + Wg{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 1);
                end
                edf = 2 * cCovBc.^2 ./ CovcCovBc - 2;
                clear CovcCovBc cCovBc
                if WB.stat == 'T'
                  if any(score > 0)
                    p(score > 0)  = spm_Tcdf(-score(score>0), edf(score>0));
                  end
                  if any(score <= 0)
                    p(score <= 0) = spm_Tcdf(score(score<=0), edf(score<=0));
                  end
                else
                  p = 2 * spm_Tcdf(-abs(score), edf);
                end
              case 3
                CovcCovBc = 0;
                for g = 1:nGr
                  CovcCovBc = CovcCovBc + Wg{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 2);
                end
                edf = 2 * cCovBc.^2 ./ CovcCovBc;
                clear CovcCovBc cCovBc
                
                if WB.stat == 'T'
                  if any(score > 0)
                    p(score > 0)  = spm_Tcdf(-score(score>0), edf(score>0));
                  end
                  if any(score <= 0)
                    p(score <= 0) = spm_Tcdf(score(score<=0), edf(score<=0));
                  end
                else
                  p = 2 * spm_Tcdf(-abs(score), edf);
                end
            end
            if SwE.WB.stat == 'F'
              score = score .^2;
            end
            activatedVoxels = [activatedVoxels, p <= WB.clusterInfo.primaryThreshold & score > 0];
            if (SwE.WB.stat == 'T')
              activatedVoxelsNeg = [activatedVoxelsNeg, p <= WB.clusterInfo.primaryThreshold & score < 0];
            end
          end
          
          maxScore(1) = max(maxScore(1), max(score));
          
          if (SwE.WB.stat == 'T')
            minScore(1) = min(minScore(1), min(score));
          end
        else
          % need to loop at every voxel
          if isfield(SwE.type,'modified')
            cCovBc = weightR * Cov_vis;
          end
          cBeta = conWB * beta;
          score = zeros(1, CrS);
          for iVox = 1:CrS
            cCovBc_vox = zeros(nSizeCon);
            cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
            cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
            score(iVox) = cBeta(:,iVox)' / cCovBc_vox * cBeta(:,iVox);
          end
          score = score / rankCon;
          maxScore(1) = max(maxScore(1), max(score));
          % save cluster information is needed
          if (SwE.WB.clusterWise == 1)
            % need to convert score into parametric p-values
            p = zeros(1, CrS);
            switch dof_type
              
              case 1
                error('degrees of freedom type still not implemented for the WB')
                
              case 2
                CovcCovBc = 0;
                for g = 1:nGr
                  CovcCovBc = CovcCovBc + Wg_testII{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 1);
                end
                edf = 2 * (sum(swe_duplication_matrix(nSizeCon), 1) * cCovBc).^2 ./ CovcCovBc - 2;
                
              case 3
                CovcCovBc = 0;
                for g = 1:nGr
                  CovcCovBc = CovcCovBc + Wg_testIII{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 2);
                end
                tmp = eye(nSizeCon);
                edf = (sum(swe_duplication_matrix(nSizeCon), 1) * cCovBc.^2 +...
                  (tmp(:)' * swe_duplication_matrix(nSizeCon) * cCovBc).^2) ./ CovcCovBc;
            end
            scoreTmp = (edf-rankCon+1) ./ edf .* score;
            scoreTmp(scoreTmp < 0 ) = 0;
            % spm_Fcdf can be inaccurate in some case --> fcdf
            if dof_type == 0
              p(scoreTmp>0) = betainc((edf-rankCon+1)./(edf-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf-rankCon+1)/2, rankCon/2);
            else
              p(scoreTmp>0) = betainc((edf(scoreTmp>0)-rankCon+1)./(edf(scoreTmp>0)-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf(scoreTmp>0)-rankCon+1)/2, rankCon/2);
              p(scoreTmp == 0) = 1;
            end
            activatedVoxels = [activatedVoxels, p <= WB.clusterInfo.primaryThreshold];
          end
        end
        
        %-Save betas etc. for current plane as we go along
        %----------------------------------------------------------
        CrYWB             = [CrYWB,    YWB]; %#ok<AGROW>
        CrResWB           = [CrResWB,  resWB]; %#ok<AGROW>
        CrScore           = [CrScore,  score]; %#ok<AGROW>
        
      end % (CrS)
      
      %-Append new inmask voxel locations and volumes
      %------------------------------------------------------------------
      XYZ(:,S + (1:CrS)) = xyz(:,Cm);     %-InMask XYZ voxel coords
      Q                  = [Q I(Cm)];     %#ok<AGROW> %-InMask XYZ voxel indices
      S                  = S + CrS;       %-Volume analysed (voxels)
      
    end % (bch)
    
    %-Plane complete, write plane to image files (unless 1st pass)
    %======================================================================
    
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...saving plane'); %-#
    
    jj = NaN(xdim,ydim,numel(CrPl));
    
    %-Write Mask image
    %------------------------------------------------------------------
    if ~isempty(Q), jj(Q) = 1; end
    VM    = spm_write_plane(VM, ~isnan(jj), CrPl);
    
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
    
    %-Write score image of the original data
    %------------------------------------------------------------------
    if ~isempty(Q), jj(Q) = CrScore; end
    Vscore = spm_write_plane(Vscore, jj, CrPl);
    
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
  clear diffVis
  
  Y     = Y(:,Cm);                          %-Data within mask
  CrS   = sum(Cm);                          %-# current voxels
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
    
    % restricted fitted data
    if WB.RWB == 1
      YWB = tmpR2 * beta;
      if SwE.WB.SS >= 4 % SC2 or SC3
        resWB = zeros(size(Y));
        for i = 1:nSubj
          resWB(iSubj==uSubj(i),:) = corrWB{i} *...
            (Y(iSubj==uSubj(i),:)-YWB(iSubj==uSubj(i),:));
        end
      else
        resWB  = diag(corrWB) * (Y-YWB);
      end
    else
      YWB = xX.X * beta;
      if SwE.WB.SS >= 4 % SC2 or SC3
        resWB = zeros(size(Y));
        for i = 1:nSubj
          resWB(iSubj==uSubj(i),:) = corrWB{i} *...
            (Y(iSubj==uSubj(i),:)-YWB(iSubj==uSubj(i),:));
        end
      else
        resWB  = diag(corrWB) * (Y-YWB);
      end
    end
    
    if WB.RSwE == 1
      if SwE.SS >= 4  % Cluster-wise adjustments
        res = zeros(size(Y));
        for i = 1:nSubj
          res(iSubj==uSubj(i),:) = corr{i} *...
            (Y(iSubj==uSubj(i),:)-tmpR2(iSubj==uSubj(i),:)*beta);
        end
      else
        res = diag(corr) * (Y - tmpR2 * beta); %-Corrected residuals
      end
    else
      if SwE.SS >= 4  % Cluster-wise adjustments
        res = zeros(size(Y));
        for i = 1:nSubj
          res(iSubj==uSubj(i),:) = corr{i} *...
            (Y(iSubj==uSubj(i),:)-xX.X(iSubj==uSubj(i),:)*beta);
        end
      else
        res = diag(corr) * (Y-xX.X*beta); %-Corrected residuals
      end
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
    else
      cCovBc = 0;
      for i = 1:nSubj
        Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
          (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
        cCovBc = cCovBc + Cov_beta_i_tmp;
      end
    end
    
    % compute the score
    if (nSizeCon == 1)
      if isfield(SwE.type,'modified')
        cCovBc = weightR * Cov_vis;
      end
      
      score = (conWB * beta) ./ sqrt(cCovBc);
      
      if (SwE.WB.clusterWise == 1)
        % need to convert score into parametric p-values
        p = zeros(1, CrS);
        switch dof_type
          case 0
            if WB.stat == 'T'
              if any(score > 0)
                p(score > 0)  = spm_Tcdf(-score(score>0), edf);
              end
              if any(score <= 0)
                p(score <= 0) = spm_Tcdf(score(score<=0), edf);
              end
            else
              p = 2 * spm_Tcdf(-abs(score), edf);
            end
          case 1
            error('degrees of freedom type still not implemented for the WB')
            
          case 2
            CovcCovBc = 0;
            for g = 1:nGr
              %CovcCovBc = CovcCovBc + Wg{g} * swe_vechCovVechV(Cov_vis(SwE.Vis.iGr_Cov_vis_g==g,:), dofMat{g}, 1);
              CovcCovBc = CovcCovBc + Wg{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 1);
            end
            edf = 2 * cCovBc.^2 ./ CovcCovBc - 2;
            clear CovcCovBc cCovBc
            if WB.stat == 'T'
              if any(score > 0)
                p(score > 0)  = spm_Tcdf(-score(score>0), edf(score>0));
              end
              if any(score <= 0)
                p(score <= 0) = spm_Tcdf(score(score<=0), edf(score<=0));
              end
            else
              p = 2 * spm_Tcdf(-abs(score), edf);
            end
          case 3
            CovcCovBc = 0;
            for g = 1:nGr
              CovcCovBc = CovcCovBc + Wg{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 2);
            end
            edf = 2 * cCovBc.^2 ./ CovcCovBc;
            clear CovcCovBc cCovBc
            
            if WB.stat == 'T'
              if any(score > 0)
                p(score > 0)  = spm_Tcdf(-score(score>0), edf(score>0));
              end
              if any(score <= 0)
                p(score <= 0) = spm_Tcdf(score(score<=0), edf(score<=0));
              end
            else
              p = 2 * spm_Tcdf(-abs(score), edf);
            end
        end
        if SwE.WB.stat == 'F'
          score = score .^2;
        end
        activatedVoxels = [activatedVoxels, p <= WB.clusterInfo.primaryThreshold & score > 0];
        if (SwE.WB.stat == 'T')
          activatedVoxelsNeg = [activatedVoxelsNeg, p <= WB.clusterInfo.primaryThreshold & score < 0];
        end
      end
      
      maxScore(1) = max(score);
      
      if (SwE.WB.stat == 'T')
        minScore(1) = min(score);
      end
    else
      % need to loop at every voxel
      if isfield(SwE.type,'modified')
        cCovBc = weightR * Cov_vis;
      end
      cBeta = conWB * beta;
      score = zeros(1, CrS);
      for iVox = 1:CrS
        cCovBc_vox = zeros(nSizeCon);
        cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
        cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
        score(iVox) = cBeta(:,iVox)' / cCovBc_vox * cBeta(:,iVox);
      end
      score = score / rankCon;
      maxScore(1) = max(score);
      % save cluster information is needed
      if (SwE.WB.clusterWise == 1)
        % need to convert score into parametric p-values
        p = zeros(1, CrS);
        switch dof_type
          
          case 1
            error('degrees of freedom type still not implemented for the WB')
            
          case 2
            CovcCovBc = 0;
            for g = 1:nGr
              CovcCovBc = CovcCovBc + Wg_testII{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 1);
            end
            edf = 2 * (sum(swe_duplication_matrix(nSizeCon), 1) * cCovBc).^2 ./ CovcCovBc - 2;
            
          case 3
            CovcCovBc = 0;
            for g = 1:nGr
              CovcCovBc = CovcCovBc + Wg_testIII{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 2);
            end
            tmp = eye(nSizeCon);
            edf = (sum(swe_duplication_matrix(nSizeCon), 1) * cCovBc.^2 +...
              (tmp(:)' * swe_duplication_matrix(nSizeCon) * cCovBc).^2) ./ CovcCovBc;
        end
        scoreTmp = (edf-rankCon+1) ./ edf .* score;
        scoreTmp(scoreTmp < 0 ) = 0;
        % spm_Fcdf can be inaccurate in some case --> fcdf
        if dof_type == 0
          p(scoreTmp>0) = betainc((edf-rankCon+1)./(edf-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf-rankCon+1)/2, rankCon/2);
        else
          p(scoreTmp>0) = betainc((edf(scoreTmp>0)-rankCon+1)./(edf(scoreTmp>0)-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf(scoreTmp>0)-rankCon+1)/2, rankCon/2);
          p(scoreTmp == 0) = 1;
        end
        activatedVoxels = [activatedVoxels, p <= WB.clusterInfo.primaryThreshold];
      end
    end
    
  end % (CrS)
  M           = [];
  DIM         = [];
  S           = CrS;  
  VM          = 'mask.mat';
  Vbeta       = 'beta.mat';
  Vscore      = 'score.mat';

  mask = Cm;       
  save('mask.mat', 'mask');
  clear mask
  
  tmp = beta;
  beta = nan(nBeta, nVox);
  beta(:,Cm) = tmp;
  save('beta.mat', 'beta');
  clear beta tmp
  
  tmp = score;
  score = nan(1, nVox);
  score(:,Cm) = tmp;
  save('score.mat', 'score');
  score = tmp;
  clear tmp
  
  
  fprintf('\n');                                                          %-#
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
if isfield(SwE.WB.clusterInfo, 'Vfaces')
  XYZ = [];
end

SwE.xVol.XYZ   = XYZ;               %-InMask XYZ coords (voxels)
SwE.xVol.M     = M;                 %-voxels -> mm
SwE.xVol.iM    = inv(M);            %-mm -> voxels
SwE.xVol.DIM   = DIM;               %-image dimensions
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
  SwE.WB.iGr_Cov_vis_g = iGr_Cov_vis_g;
end
SwE.VM         = VM;                %-Filehandle - Mask

SwE.xX         = xX;                %-design structure
SwE.xM         = xM;                %-mask structure

SwE.swd        = pwd;

SwE.Subj.uSubj = uSubj;
SwE.Subj.nSubj = nSubj;

if isfield(SwE.type,'modified')
  
  SwE.Vis.uVis_g = uVis_g;
  SwE.Vis.nVis_g = nVis_g;
  SwE.Vis.nCov_vis_g = nCov_vis_g;
  SwE.Vis.nCov_vis = nCov_vis;
  
  SwE.Gr.uGr       = uGr;
  SwE.Gr.nGr       = nGr;
  SwE.Gr.nSubj_g   = nSubj_g;
  SwE.Gr.uSubj_g   = uSubj_g;
  
end

%-Save analysis parameters in SwE.mat file
%--------------------------------------------------------------------------
if spm_matlab_version_chk('7') >=0
  save('SwE','SwE','-V6');
else
  save('SwE','SwE');
end

%-Save analysis original max min in files
%--------------------------------------------------------------------------
save('maxScore.mat', 'maxScore');
if (SwE.WB.clusterWise == 1)
  save('maxClusterSize.mat', 'maxClusterSize');
end
if (SwE.WB.stat == 'T')
  save('minScore.mat', 'minScore');
  if (SwE.WB.clusterWise == 1)
    save('maxClusterSizeNeg.mat', 'maxClusterSizeNeg');
  end
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
      if WB.RSwE == 0 % U-SwE
        if SwE.SS >= 4  % Cluster-wise adjustments SC2 or SC3
          res = zeros(size(Y_b));
          for i = 1:nSubj
            res(iSubj==uSubj(i),:) = corr{i} *...
              (Y_b(iSubj==uSubj(i),:) - xX.X(iSubj==uSubj(i),:) * beta);
          end
        else
          res = diag(corr) * (Y_b - xX.X * beta); %-Corrected residuals
        end
      else % R-SwE
        if SwE.SS >= 4  % Cluster-wise adjustments SC2 or SC3
          res = zeros(size(Y_b));
          for i = 1:nSubj
            res(iSubj==uSubj(i),:) = corr{i} *...
              (Y_b(iSubj==uSubj(i),:) - tmpR2(iSubj==uSubj(i),:) * beta);
          end
        else
          res = diag(corr) * (Y_b - tmpR2 * beta); %-Corrected residuals
        end
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
      else
        cCovBc = 0;
        for i = 1:nSubj
          Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
            (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
          cCovBc = cCovBc + Cov_beta_i_tmp;
        end
      end
      
      % compute the score
      if (nSizeCon == 1)
        if isfield(SwE.type,'modified')
          cCovBc = weightR * Cov_vis;
        end
        score = (conWB * beta) ./ sqrt(cCovBc);
        
        clear beta
        
        if (WB.clusterWise == 1)
          % need to convert score into parametric p-values
          p = zeros(1, blksz);
          switch dof_type
            case 0
              if WB.stat == 'T'
                if any(score > 0)
                  p(score > 0)  = spm_Tcdf(-score(score>0), edf);
                end
                if any(score <= 0)
                  p(score <= 0) = spm_Tcdf(score(score<=0), edf);
                end
              else
                p = 2 * spm_Tcdf(-abs(score), edf);
              end
            case 1
              error('degrees of freedom type still not implemented for the WB')
              
            case 2
              CovcCovBc = 0;
              for g = 1:nGr
                CovcCovBc = CovcCovBc + Wg{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 1);
              end
              edf = 2 * cCovBc.^2 ./ CovcCovBc - 2;
              clear CovcCovBc cCovBc
              
              if WB.stat == 'T'
                if any(score > 0)
                  p(score > 0)  = spm_Tcdf(-score(score>0), edf(score>0));
                end
                if any(score <= 0)
                  p(score <= 0) = spm_Tcdf(score(score<=0), edf(score<=0));
                end
              else
                p = 2 * spm_Tcdf(-abs(score), edf);
              end
            case 3
              CovcCovBc = 0;
              for g = 1:nGr
                CovcCovBc = CovcCovBc + Wg{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 2);
              end
              edf = 2 * cCovBc.^2 ./ CovcCovBc;
              clear CovcCovBc cCovBc
              
              if WB.stat == 'T'
                if any(score > 0)
                  p(score > 0)  = spm_Tcdf(-score(score>0), edf(score>0));
                end
                if any(score <= 0)
                  p(score <= 0) = spm_Tcdf(score(score<=0), edf(score<=0));
                end
              else
                p = 2 * spm_Tcdf(-abs(score), edf);
              end
          end
          if WB.stat == 'F'
            score = score .^2;
          end
          activatedVoxels(index) = p <= WB.clusterInfo.primaryThreshold & score > 0;
          if (WB.stat == 'T')
            activatedVoxelsNeg(index) = p <= WB.clusterInfo.primaryThreshold & score < 0;
          end
        end
        
        
        uncP(index) = uncP(index) + (score >= originalScore(index)) * 1;
        
        maxScore(b+1) = max(maxScore(b+1), max(score));
        if (WB.stat == 'T')
          minScore(b+1) = min(minScore(b+1), min(score));
        end
      else
        % need to loop at every voxel
        if isfield(SwE.type,'modified')
          cCovBc = weightR * Cov_vis;
        end
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
        % save cluster information is needed
        if (WB.clusterWise == 1)
          % need to convert score into parametric p-values
          p = zeros(1, blksz);
          switch dof_type
            
            case 1
              error('degrees of freedom type still not implemented for the WB')
              
            case 2
              CovcCovBc = 0;
              for g = 1:nGr
                CovcCovBc = CovcCovBc + Wg_testII{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 1);
              end
              edf = 2 * (sum(swe_duplication_matrix(nSizeCon), 1) * cCovBc).^2 ./ CovcCovBc - 2;
              
            case 3
              CovcCovBc = 0;
              for g = 1:nGr
                CovcCovBc = CovcCovBc + Wg_testIII{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 2);
              end
              tmp = eye(nSizeCon);
              edf = (sum(swe_duplication_matrix(nSizeCon), 1) * cCovBc.^2 +...
                (tmp(:)' * swe_duplication_matrix(nSizeCon) * cCovBc).^2) ./ CovcCovBc;
              
          end
          scoreTmp = (edf-rankCon+1) ./ edf .* score;
          scoreTmp(scoreTmp < 0 ) = 0;
          if dof_type == 0
            p(scoreTmp>0) = betainc((edf-rankCon+1)./(edf-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf-rankCon+1)/2, rankCon/2);
          else
            p(scoreTmp>0) = betainc((edf(scoreTmp>0)-rankCon+1)./(edf(scoreTmp>0)-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf(scoreTmp>0)-rankCon+1)/2, rankCon/2);
            p(scoreTmp == 0) = 1;
          end
          
          activatedVoxels(index) = p <= WB.clusterInfo.primaryThreshold;
        end
        uncP(index) = uncP(index) + (score >= originalScore(index)) * 1;
        
        maxScore(b+1) = max(maxScore(b+1), max(score));
      end
      
    end % (bch)
  else
    
    %-Print progress information in command window
    %------------------------------------------------------------------
    str = sprintf('Bootstrap # %i', b);
    
    fprintf('%-40s: %1s',str,' ');
    
    Y_b = YWB + resWB .* repmat(resamplingMatrix(:,b),1,S);
    
    beta  = pX * Y_b;                     %-Parameter estimates
    if WB.RSwE == 0 % U-SwE
      if SwE.SS >= 4  % Cluster-wise adjustments SC2 or SC3
        res = zeros(size(Y_b));
        for i = 1:nSubj
          res(iSubj==uSubj(i),:) = corr{i} *...
            (Y_b(iSubj==uSubj(i),:) - xX.X(iSubj==uSubj(i),:) * beta);
        end
      else
        res = diag(corr) * (Y_b - xX.X * beta); %-Corrected residuals
      end
    else % R-SwE
      if SwE.SS >= 4  % Cluster-wise adjustments SC2 or SC3
        res = zeros(size(Y_b));
        for i = 1:nSubj
          res(iSubj==uSubj(i),:) = corr{i} *...
            (Y_b(iSubj==uSubj(i),:) - tmpR2(iSubj==uSubj(i),:) * beta);
        end
      else
        res = diag(corr) * (Y_b - tmpR2 * beta); %-Corrected residuals
      end
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
    else
      cCovBc = 0;
      for i = 1:nSubj
        Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
          (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
        cCovBc = cCovBc + Cov_beta_i_tmp;
      end
    end
    
    % compute the score
    if (nSizeCon == 1)
      if isfield(SwE.type,'modified')
        cCovBc = weightR * Cov_vis;
      end
      score = (conWB * beta) ./ sqrt(cCovBc);
      
      clear beta
      
      if (WB.clusterWise == 1)
        % need to convert score into parametric p-values
        p = zeros(1, S);
        switch dof_type
          case 0
            if WB.stat == 'T'
              if any(score > 0)
                p(score > 0)  = spm_Tcdf(-score(score>0), edf);
              end
              if any(score <= 0)
                p(score <= 0) = spm_Tcdf(score(score<=0), edf);
              end
            else
              p = 2 * spm_Tcdf(-abs(score), edf);
            end
          case 1
            error('degrees of freedom type still not implemented for the WB')
            
          case 2
            CovcCovBc = 0;
            for g = 1:nGr
              CovcCovBc = CovcCovBc + Wg{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 1);
            end
            edf = 2 * cCovBc.^2 ./ CovcCovBc - 2;
            clear CovcCovBc cCovBc
            
            if WB.stat == 'T'
              if any(score > 0)
                p(score > 0)  = spm_Tcdf(-score(score>0), edf(score>0));
              end
              if any(score <= 0)
                p(score <= 0) = spm_Tcdf(score(score<=0), edf(score<=0));
              end
            else
              p = 2 * spm_Tcdf(-abs(score), edf);
            end
          case 3
            CovcCovBc = 0;
            for g = 1:nGr
              CovcCovBc = CovcCovBc + Wg{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 2);
            end
            edf = 2 * cCovBc.^2 ./ CovcCovBc;
            clear CovcCovBc cCovBc
            
            if WB.stat == 'T'
              if any(score > 0)
                p(score > 0)  = spm_Tcdf(-score(score>0), edf(score>0));
              end
              if any(score <= 0)
                p(score <= 0) = spm_Tcdf(score(score<=0), edf(score<=0));
              end
            else
              p = 2 * spm_Tcdf(-abs(score), edf);
            end
        end
        if WB.stat == 'F'
          score = score .^2;
        end
        activatedVoxels = p <= WB.clusterInfo.primaryThreshold & score > 0;
        if (WB.stat == 'T')
          activatedVoxelsNeg = p <= WB.clusterInfo.primaryThreshold & score < 0;
        end
      end
      
      uncP = uncP + (score >= originalScore) * 1;
      
      maxScore(b+1) = max(score);
      if (WB.stat == 'T')
        minScore(b+1) = min(score);
      end
    else
      % need to loop at every voxel
      if isfield(SwE.type,'modified')
        cCovBc = weightR * Cov_vis;
      end
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
      % save cluster information is needed
      if (WB.clusterWise == 1)
        % need to convert score into parametric p-values
        p = zeros(1, S);
        switch dof_type
          
          case 1
            error('degrees of freedom type still not implemented for the WB')
            
          case 2
            CovcCovBc = 0;
            for g = 1:nGr
              CovcCovBc = CovcCovBc + Wg_testII{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 1);
            end
            edf = 2 * (sum(swe_duplication_matrix(nSizeCon), 1) * cCovBc).^2 ./ CovcCovBc - 2;
            
          case 3
            CovcCovBc = 0;
            for g = 1:nGr
              CovcCovBc = CovcCovBc + Wg_testIII{g} * swe_vechCovVechV(Cov_vis(iGr_Cov_vis_g==g,:), dofMat{g}, 2);
            end
            tmp = eye(nSizeCon);
            edf = (sum(swe_duplication_matrix(nSizeCon), 1) * cCovBc.^2 +...
              (tmp(:)' * swe_duplication_matrix(nSizeCon) * cCovBc).^2) ./ CovcCovBc;
            
        end
        scoreTmp = (edf-rankCon+1) ./ edf .* score;
        scoreTmp(scoreTmp < 0 ) = 0;
        if dof_type == 0
          p(scoreTmp>0) = betainc((edf-rankCon+1)./(edf-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf-rankCon+1)/2, rankCon/2);
        else
          p(scoreTmp>0) = betainc((edf(scoreTmp>0)-rankCon+1)./(edf(scoreTmp>0)-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf(scoreTmp>0)-rankCon+1)/2, rankCon/2);
          p(scoreTmp == 0) = 1;
        end
        
        activatedVoxels = p <= WB.clusterInfo.primaryThreshold;
      end
      uncP = uncP + (score >= originalScore) * 1;
      
      maxScore(b+1) = max(score);
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
  %-Save analysis original max min in files
  %--------------------------------------------------------------------------
  save('maxScore.mat', 'maxScore');
  if (WB.clusterWise == 1)
    save('maxClusterSize.mat', 'maxClusterSize');
  end
  if (WB.stat == 'T')
    save('minScore.mat', 'minScore');
    if (WB.clusterWise == 1)
      save('maxClusterSizeNeg.mat', 'maxClusterSizeNeg');
    end
  end
  toc
  spm_progress_bar('Set',100 * b / WB.nB);
end


%==========================================================================
%- produce results images
%==========================================================================
if isMat
  uncP = uncP / (WB.nB + 1);
  uncP_pos = nan(1, nVox);
  uncP_pos(:,Cm) = uncP;
  save('uncP_pos.mat', 'uncP_pos');
  lUncP_pos = -log10(uncP);
  save('lUncP_pos.mat', 'lUncP_pos');
  clear lUncP_pos
  
  
  if WB.stat == 'T'
    uncP_neg = 1 + 1/(WB.nB + 1) - uncP_pos;
    save('uncP_neg.mat', 'uncP_neg');
    lUncP_neg = -log10(uncP_neg);
    save('lUncP_neg.mat', 'lUncP_neg');
    clear lUncP_neg
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
    FWERP = FWERP + (maxScore(b+1) > originalScore - tol) * 1;
  end
  FWERP = FWERP / (WB.nB + 1);
  fwerP_pos = nan(1, nVox);
  fwerP_pos(:,Cm) = FWERP;
  save('fwerP_pos.mat', 'fwerP_pos');
  lFwerP_pos = -log10(fwerP_pos);
  save('lFwerP_pos.mat', 'lFwerP_pos');
  clear lFwerP_pos fwerP_pos FWERP
  
  
  if WB.stat == 'T'
    FWERPNeg = ones(1, S); % 1 because the original maxScore is always > original Score
    
    for b = 1:WB.nB
      %-FWER-corrected p is proportion of randomisation greater or
      % equal to statistic.
      %-Use a > b -tol rather than a >= b to avoid comparing
      % two reals for equality.
      FWERPNeg = FWERPNeg + (minScore(b+1) < originalScore + tol) * 1;
    end
    FWERPNeg = FWERPNeg / (WB.nB + 1);
    fwerP_neg = nan(1, nVox);
    fwerP_neg(:,Cm) = FWERPNeg;
    save('fwerP_neg.mat', 'fwerP_neg');
    lFwerP_neg = -log10(fwerP_neg);
    save('lFwerP_neg.mat', 'lFwerP_neg');
    clear lFwerP_neg fwerP_neg
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
  save('fdrP_pos.mat', 'fdrP_pos');
  lFdrP_pos = -log10(fdrP_pos);
  save('lFdrP_pos.mat', 'lFdrP_pos');
  clear lFdrP_pos fdrP_pos fdrP
  
  if WB.stat =='T'
    try
      fdrP = spm_P_FDR(1 + 1/(WB.nB + 1) - uncP);
    catch
      fdrP = spm_P_FDR(1 + 1/(WB.nB + 1) - uncP,[],'P',[],sort(1 + 1/(WB.nB + 1) - uncP)');
    end
    fdrP_neg = nan(1, nVox);
    fdrP_neg(:,Cm) = fdrP;
    save('fdrP_neg.mat', 'fdrP_neg');
    lFdrP_neg = -log10(fdrP_neg);
    save('lFdrP_neg.mat', 'lFdrP_neg');
    clear lFdrP_neg fdrP_neg fdrP
  end
  
  if WB.clusterWise == 1
    % Not sure what to output. So might be changed later.
    % For now, -log(p_{cluster-wise FWER}) image with nan for non-surviving
    % voxels after the thresholding of the original data
    
    clusterFWERP = ones(1, SwE.WB.clusterInfo.nCluster); % 1 because the original maxScore is always > original Score
    if (~isempty(SwE.WB.clusterInfo.clusterSize))
      for b = 1:WB.nB
        clusterFWERP = clusterFWERP + (maxClusterSize(b+1) >= SwE.WB.clusterInfo.clusterSize) * 1;
      end
      clusterFWERP = clusterFWERP / (WB.nB + 1);
    end
    clusterFwerP_pos_perCluster   = clusterFWERP;
    lClusterFwerP_pos_perCluster  = -log10(clusterFwerP_pos_perCluster);
    save('clusterFwerP_pos_perCluster.mat', 'clusterFwerP_pos_perCluster');
    save('lClusterFwerP_pos_perCluster.mat', 'lClusterFwerP_pos_perCluster');
    
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
    lClusterFwerP_pos_perElement  = -log10(clusterFwerP_pos_perElement);
    save('clusterFwerP_pos_perElement.mat', 'clusterFwerP_pos_perElement');
    save('lClusterFwerP_pos_perElement.mat', 'lClusterFwerP_pos_perElement');
    
    if WB.stat =='T'
      
      clusterFWERPNeg = ones(1, SwE.WB.clusterInfo.nClusterNeg); % 1 because the original maxScore is always > original Score
      if (~isempty(SwE.WB.clusterInfo.clusterSizeNeg))
        for b = 1:WB.nB
          clusterFWERPNeg = clusterFWERPNeg + (maxClusterSizeNeg(b+1) >= SwE.WB.clusterInfo.clusterSizeNeg) * 1;
        end
        clusterFWERPNeg = clusterFWERPNeg / (WB.nB + 1);
      end
      clusterFwerP_neg_perCluster   = clusterFWERPNeg;
      lClusterFwerP_neg_perCluster  = -log10(clusterFwerP_neg_perCluster);
      save('clusterFwerP_neg_perCluster.mat', 'clusterFwerP_neg_perCluster');
      save('lClusterFwerP_neg_perCluster.mat', 'lClusterFwerP_neg_perCluster');
      
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
      lClusterFwerP_neg_perElement  = -log10(clusterFwerP_neg_perElement);
      save('clusterFwerP_neg_perElement.mat', 'clusterFwerP_neg_perElement');
      save('lClusterFwerP_neg_perElement.mat', 'lClusterFwerP_neg_perElement');
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
  spm_write_vol(VlP_pos, tmp);
  
  if WB.stat == 'T'
    tmp(Q) = -log10(1 + 1/(WB.nB + 1) - uncP);
    spm_write_vol(VlP_neg, tmp);
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
    FWERP = FWERP + (maxScore(b+1) > originalScore - tol) * 1;
  end
  FWERP = FWERP / (WB.nB + 1);
  tmp(Q) = -log10(FWERP);
  spm_write_vol(VlP_FWE_pos, tmp);
  
  if WB.stat == 'T'
    FWERPNeg = ones(1, S); % 1 because the original maxScore is always > original Score
    
    for b = 1:WB.nB
      %-FWER-corrected p is proportion of randomisation greater or
      % equal to statistic.
      %-Use a > b -tol rather than a >= b to avoid comparing
      % two reals for equality.
      FWERPNeg = FWERPNeg + (minScore(b+1) < originalScore + tol) * 1;
    end
    FWERPNeg = FWERPNeg / (WB.nB + 1);
    tmp(Q) = -log10(FWERPNeg);
    spm_write_vol(VlP_FWE_neg, tmp);
  end
  
  %
  % - write out lP_FDR+ and lP_FDR- images;
  %
  try
    tmp(Q) = -log10(spm_P_FDR(uncP));
  catch
    tmp(Q) = -log10(spm_P_FDR(uncP,[],'P',[],sort(uncP)'));
  end
  spm_write_vol(VlP_FDR_pos, tmp);
  
  if WB.stat =='T'
    try
      tmp(Q) = -log10(spm_P_FDR(1 + 1/(WB.nB + 1) - uncP));
    catch
      tmp(Q) = -log10(spm_P_FDR(1 + 1/(WB.nB + 1) - uncP,[],'P',[],sort(1 + 1/(WB.nB + 1) - uncP)'));
    end
    spm_write_vol(VlP_FDR_neg, tmp);
  end
  
  if WB.clusterWise == 1
    % Not sure what to output. So might be changed later.
    % For now, -log(p_{cluster-wise FWER}) image with nan for non-surviving
    % voxels after the thresholding of the original data
    Q = cumprod([1,SwE.xVol.DIM(1:2)']) * SwE.WB.clusterInfo.LocActivatedVoxels - ...
      sum(cumprod(SwE.xVol.DIM(1:2)'));
    tmp= nan(SwE.xVol.DIM');
    
    clusterFWERP = ones(1, SwE.WB.clusterInfo.nCluster); % 1 because the original maxScore is always > original Score
    if (~isempty(SwE.WB.clusterInfo.clusterSize))
      for b = 1:WB.nB
        clusterFWERP = clusterFWERP + (maxClusterSize(b+1) >= SwE.WB.clusterInfo.clusterSize) * 1;
      end
      clusterFWERP = clusterFWERP / (WB.nB + 1);
    end
    tmp2 = -log10(clusterFWERP);
    
    tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxels,2));
    for iC = 1:SwE.WB.clusterInfo.nCluster
      tmp3(SwE.WB.clusterInfo.clusterAssignment == iC) = tmp2(iC);
    end
    tmp(Q) = tmp3;
    spm_write_vol(VlP_clusterFWE_pos, tmp);
    if WB.stat =='T'
      Q = cumprod([1,SwE.xVol.DIM(1:2)']) * SwE.WB.clusterInfo.LocActivatedVoxelsNeg - ...
        sum(cumprod(SwE.xVol.DIM(1:2)'));
      tmp= nan(SwE.xVol.DIM');
      
      clusterFWERPNeg = ones(1, SwE.WB.clusterInfo.nClusterNeg); % 1 because the original maxScore is always > original Score
      if (~isempty(SwE.WB.clusterInfo.clusterSizeNeg))
        for b = 1:WB.nB
          clusterFWERPNeg = clusterFWERPNeg + (maxClusterSizeNeg(b+1) >= SwE.WB.clusterInfo.clusterSizeNeg) * 1;
        end
        clusterFWERPNeg = clusterFWERPNeg / (WB.nB + 1);
      end
      tmp2 = -log10(clusterFWERPNeg);
      
      tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxelsNeg, 2));
      for iC = 1:SwE.WB.clusterInfo.nClusterNeg
        tmp3(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = tmp2(iC);
      end
      tmp(Q) = tmp3;
      spm_write_vol(VlP_clusterFWE_neg, tmp);
    end
  end
end
%==========================================================================
%- E N D: Cleanup GUI
%==========================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#
%spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
fprintf('...use the saved images for assessment\n\n')
