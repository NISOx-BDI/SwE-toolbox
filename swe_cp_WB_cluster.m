function swe_cp_WB_cluster(SwE)

%-Say hello
%--------------------------------------------------------------------------
Finter = spm('CreateIntWin','off');
set(Finter,'name','SwE estimation');
set(Finter,'vis','on')

%-Get SwE.mat[s] if necessary
%--------------------------------------------------------------------------
if nargin == 0
  P     = cellstr(spm_select(Inf,'^SwE\.mat$','Select SwE.mat[s]'));
  for i = 1:length(P)
    swd     = fileparts(P{i});
    load(fullfile(swd,'SwE.mat'));
    SwE.swd = swd;
    swe_cp_WB_cluster(SwE);
  end
  return
end
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

%-Delete files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(SwE.swd,'mask.img'),'file') == 2
  
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

files = {'^mask\..{3}$','^beta_.{4}\..{3}$','^con_.{4}\..{3}$',...
  '^ResI_.{4}\..{3}$','^cov_beta_.{4}d_.{4}\..{3}$',...
  '^cov_beta_g_.{4}d_.{4}d_.{4}\..{3}$',...
  '^cov_vis_.{4}d_.{4}d_.{4}\..{3}$'
  };

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
conWB = SwE.WB.con;
nSizeCon = size(conWB,1);
rankCon = rank(conWB);
nCov_betaR = (nSizeCon + 1) * nSizeCon / 2;
tmpR = (xX.X' * xX.X) \ conWB';
tmpR = tmpR / (conWB * tmpR);
tmpR2 = xX.X * (eye(nBeta) - tmpR * conWB);
HatR = xX.X * (pX - tmpR * conWB * pX);

% only for the HCP data
edf = 79;

% score multiplicator value if contrast of rank > 1
if (nSizeCon > 1)
  scoreMult = (edf-rankCon+1)/edf/rankCon;
end

switch SwE.WB.SS
  case 2
    corrR  = (1-diag(HatR)).^(-0.5); % residual correction (type 2)
  case 3
    corrR  = (1-diag(HatR)).^(-1); % residual correction (type 3)
  case 4
    corrR =cell(nSubj,1);
    I_Hat = eye(nScan) - HatR;
    for i = 1:nSubj
      tmp = I_Hat(iSubj==uSubj(i), iSubj==uSubj(i));
      tmp = (tmp + tmp')/2;
      [tmpV, tmpE] = eig(tmp);
      corrR{i} = tmpV * diag(1./sqrt(diag(tmpE))) * tmpV';
    end
    clear I_Hat tmp
end

if (SwE.WB.RSwE == 0)
  Hat = xX.X*(pX); % Hat matrix 
  switch SwE.SS   
      case 1
          corr  = sqrt(nScan/(nScan-nBeta)); % residual correction (type 1) 
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
  for i = find(iGr_dof==0)
    if any(xX.X(i,:) & xX.X(scan,:))
      iGr_dof(i) = it;
    end
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
  if nGr_dof==1 | all(sum(tmp)==1) %#ok<OR2>
    break % all is ok, just stop the while
  else
    ind1 = find(sum(tmp)>1,1); % detect the first column in common
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

%-preprocessing for the modified SwE
if isfield(SwE.type,'modified')
  iVis      = SwE.Vis.iVis;
  iGr       = SwE.Gr.iGr;
  uGr       = unique(iGr);
  nGr       = length(uGr);
  
  % info specific for each group
  uVis_g = cell(1,nGr); % unique visits for each group
  nVis_g = zeros(1,nGr); % number of visits for each group
  uSubj_g = cell(1,nGr); % unique visits for each group
  nSubj_g = zeros(1,nGr); % number of visits for each group
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
YNaNrep  = 0

%-Maximum number of residual images for smoothness estimation
%--------------------------------------------------------------------------
MAXRES   = Inf;
nSres    = nScan;

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');               %-#

fprintf('%-40s: %30s','Output images','...initialising');           %-#

%-Initialise new mask name: current mask & conditions on voxels
%----------------------------------------------------------------------
VM    = struct('fname',  'mask.img',...
  'dim',    DIM',...
  'dt',     [spm_type('uint8') spm_platform('bigend')],...
  'mat',    M,...
  'pinfo',  [1 0 0]',...
  'descrip','swe_cp_WB:resultant analysis mask');
VM    = spm_create_vol(VM);

%-Initialise beta image files
%----------------------------------------------------------------------

Vscore = deal(struct(...
  'fname',    [],...
  'dim',      DIM',...
  'dt',       [spm_type('float32') spm_platform('bigend')],...
  'mat',      M,...
  'pinfo',    [1 0 0]',...
  'descrip',  ''));

Vscore.fname   = sprintf('score.img');
Vscore.descrip = sprintf('score original data');
Vscore = spm_create_vol(Vscore);

%-Initialise standardised residual images
%----------------------------------------------------------------------
VResR(1:nScan) = deal(struct(...
  'fname',    [],...
  'dim',      DIM',...
  'dt',       [spm_type('float32') spm_platform('bigend')],...
  'mat',      M,...
  'pinfo',    [1 0 0]',...
  'descrip',  'swe_cp_WB:StandardisedResiduals'));

for i = 1:nScan
  VResR(i).fname   = sprintf('ResR_%04d.img', i);
  VResR(i).descrip = sprintf('adjusted restricted residuals (%04d)', i);
end
VResR = spm_create_vol(VResR);
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');    %-#

%-Initialise standardised residual images
%----------------------------------------------------------------------
VYR(1:nScan) = deal(struct(...
  'fname',    [],...
  'dim',      DIM',...
  'dt',       [spm_type('float32') spm_platform('bigend')],...
  'mat',      M,...
  'pinfo',    [1 0 0]',...
  'descrip',  'swe_cp_WB:StandardisedResiduals'));

for i = 1:nScan
  VYR(i).fname   = sprintf('YfittedR_%04d.img', i);
  VYR(i).descrip = sprintf('fitted data under H0 (%04d)', i);
end
VYR = spm_create_vol(VYR);
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
if (SwE.WB.clusterWise == 1)
  activatedVoxels = false(0);
  if (SwE.WB.stat == 'T')
    activatedVoxelsNeg = false(0);
  end
end
maxScore = nan;
minScore = nan;

for z = 1:nbz:zdim                       %-loop over planes (2D or 3D data)
  
  % current plane-specific parameters
  %----------------------------------------------------------------------
  CrPl         = z:min(z+nbz-1,zdim);       %-plane list
  zords        = CrPl(:)*ones(1,xdim*ydim); %-plane Z coordinates
  CrScore      = [];                        %-scores
  CrYR         = [];                        %-fitted data under H0
  CrResR       = [];                        %-residuals
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
    % matrix design
    %------------------------------------------------------------------
    for g = 1:nGr_dof
      Cm(Cm) = any(diff(Y(iGr_dof==g,Cm),1));
    end
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
      YR = tmpR2 * beta;

      if SwE.SS == 4 % SC2
        resR = zeros(size(Y));
        for i = 1:nSubj
          resR(iSubj==uSubj(i),:) = corrR{i} *...
            (Y(iSubj==uSubj(i),:)-YR(iSubj==uSubj(i),:));
        end
        % unrestricted residuals if needed
        if (SwE.WB.RSwE == 0)
          res = zeros(size(Y));
          for i = 1:nSubj
          	res(iSubj==uSubj(i),:) = corr{i} *...
            	(Y(iSubj==uSubj(i),:)-xX.X(iSubj==uSubj(i),:)*beta);
          end 
        end
      else
        resR  = (Y-YR) .* repmat(corrR, 1, CrS);
        if (SwE.WB.RSwE == 0)
          res = diag(corr)*(Y-xX.X*beta); %-Corrected residuals
        end
      end
      clear Y                           %-Clear to save memory
      %-Estimation of the data variance-covariance components (modified SwE)
      %-SwE estimation (classic version)
      %--------------------------------------------------------------
      if isfield(SwE.type,'modified')
        Cov_vis=zeros(nCov_vis,CrS);
        if (SwE.WB.RSwE == 1) % R-SwE
          for i = Ind_Cov_vis_diag
            Cov_vis(i,:) = mean(resR(Flagk(i,:),:).^2);
          end
        else % U-SwE
          for i = Ind_Cov_vis_diag
            Cov_vis(i,:) = mean(res(Flagk(i,:),:).^2);
          end
        end
        % Check if some voxels have variance = 0 and mask them
        tmp = ~any(Cov_vis(Ind_Cov_vis_diag,:)==0);
        if ~tmp
          beta    = beta(tmp);
          resR    = resR(tmp);
          res     = res(tmp);
          YR      = YR(tmp);
          Cm(Cm)  = tmp;
          CrS     = sum(Cm);
          Cov_vis = Cov_vis(tmp);
        end
        if CrS % Check if there is at least one voxel left
          if (SwE.WB.RSwE == 1) % R-SwE
            for i = Ind_Cov_vis_off_diag
              Cov_vis(i,:)= sum(resR(Flagk(i,:),:).*resR(Flagkk(i,:),:)).*...
                sqrt(Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,1)),:).*...
                Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,2)),:)./...
                sum(resR(Flagk(i,:),:).^2)./...
                sum(resR(Flagkk(i,:),:).^2));
            end
          else % U-SwE
            for i = Ind_Cov_vis_off_diag
              Cov_vis(i,:)= sum(res(Flagk(i,:),:).*res(Flagkk(i,:),:)).*...
                sqrt(Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,1)),:).*...
                Cov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,2)),:)./...
                sum(res(Flagk(i,:),:).^2)./...
                sum(res(Flagkk(i,:),:).^2));
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
      end
      if (nSizeCon == 1)
        if (SwE.WB.stat == 'T')
          score = (conWB * beta) ./ sqrt(weightR * Cov_vis);
          minScore = min(minScore, min(score));
        else
          score = (conWB * beta).^2 ./ (weightR * Cov_vis);
        end
        maxScore = max(maxScore, max(score));
        
        if (SwE.WB.clusterWise == 1)
          activatedVoxels = [activatedVoxels, score >= SwE.WB.primaryThreshold];
          if (SwE.WB.stat == 'T')
            activatedVoxelsNeg = [activatedVoxelsNeg, score <= -SwE.WB.primaryThreshold];
          end
        end
      else
        % need to loop at every voxel
        cCovBc = weightR * Cov_vis;
        cBeta = conWB * beta;
        score = zeros(1, CrS);
        for iVox = 1:CrS
          cCovBc_vox = zeros(nSizeCon);
          cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
          cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
          score(iVox) = cBeta(:,iVox)' / cCovBc_vox * cBeta(:,iVox);
          score(iVox) = scoreMult * score(iVox); % need to multiply the score as described by the test
        end
        % save cluster information is needed
        if (SwE.WB.clusterWise == 1)
          activatedVoxels = [activatedVoxels, score >= SwE.WB.primaryThreshold];
        end
        maxScore = max(maxScore, max(score));
      end
      
      %-Save betas etc. for current plane as we go along
      %----------------------------------------------------------
      CrYR              = [CrYR,    YR]; %#ok<AGROW>
      CrResR            = [CrResR,  resR]; %#ok<AGROW>
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
  
  %-Write fitted data under H0 images
  %------------------------------------------------------------------
  for i = 1:nScan
    if ~isempty(Q), jj(Q) = CrYR(i,:); end
    VYR(i) = spm_write_plane(VYR(i), jj, CrPl);
  end
  
  %-Write fitted data under H0 images
  %------------------------------------------------------------------
  for i = 1:nScan
    if ~isempty(Q), jj(Q) = CrResR(i,:); end
    VResR(i) = spm_write_plane(VResR(i), jj, CrPl);
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
clear beta res Cov_vis CrYR CrResR CrCov_vis jj%-Clear to save memory

XYZ   = XYZ(:,1:S); % remove all the data not used

% compute the max cluster size if needed (so many ways this can be
% done... Not sure this solution is the best)
if (SwE.WB.clusterWise == 1)
  LocActivatedVoxels = XYZ(:,activatedVoxels);
  clusterAssignment = spm_clusters(LocActivatedVoxels);
  nCluster     = max(clusterAssignment);
  clusterSize = histc(clusterAssignment,1:nCluster);
  maxClusterSize = max(clusterSize);
  if (SwE.WB.stat == 'T')
    LocActivatedVoxelsNeg = XYZ(:,activatedVoxelsNeg);
    clusterAssignmentNeg = spm_clusters(LocActivatedVoxelsNeg);
    nClusterNeg     = max(clusterAssignmentNeg);
    clusterSizeNeg = histc(clusterAssignmentNeg,1:nClusterNeg);
    maxClusterSizeNeg = max(clusterSizeNeg);    
  end
end
%==========================================================================
% - P O S T   E S T I M A T I O N   C L E A N U P
%==========================================================================
if S == 0, spm('alert!','No inmask voxels - empty analysis!'); return; end


%-Smoothness estimates of component fields and RESEL counts for volume
%==========================================================================
% try
%     FWHM = SwE.xVol.FWHM;
%     VRpv = SwE.xVol.VRpv;
%     R    = SwE.xVol.R;
% catch
%     erdf      = spm_SpUtil('trRV',xX.X); % Working error df / do not agree to be checked
%     [FWHM,VRpv,R] = spm_est_smoothness(VResI,VM,[nScan erdf]);
% end

%-Compute scaled design matrix for display purposes
%--------------------------------------------------------------------------
%xX.nX        = spm_DesMtx('sca',xX,xX.name);

%-Save remaining results files and analysis parameters
%==========================================================================
fprintf('%-40s: %30s','Saving results','...writing');

%-place fields in SwE
%--------------------------------------------------------------------------
SwE.xVol.XYZ   = XYZ;               %-InMask XYZ coords (voxels)
SwE.xVol.M     = M;                 %-voxels -> mm
SwE.xVol.iM    = inv(M);            %-mm -> voxels
SwE.xVol.DIM   = DIM;               %-image dimensions
SwE.xVol.units = {'mm' 'mm' 'mm'};

SwE.WB.VYR        = VYR;               %-Filehandle - fitted data under H0
SwE.WB.VResR      = VResR;             %-Filehandle - adjusted resticted residuals
SwE.WB.Vscore     = Vscore;            %-Filehandle - score original image

SwE.WB.weightR    = weightR;
if (SwE.WB.RSwE == 1)
  SwE.WB.corrR      = corrR;
  SwE.WB.tmpR2      = tmpR2;
else
  SwE.WB.corr      = corr;
end

% cluster-wise specific fields if needed
if (SwE.WB.clusterWise == 1)
  SwE.WB.clusterInfo = [];
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
SwE.WB.Ind_Cov_vis_diag = Ind_Cov_vis_diag;
SwE.WB.Ind_Cov_vis_off_diag = Ind_Cov_vis_off_diag;
SwE.WB.Ind_corr_diag = Ind_corr_diag;
SwE.WB.Flagk      = Flagk;
SwE.WB.Flagkk     = Flagkk;
SwE.WB.iGr_Cov_vis_g = iGr_Cov_vis_g;

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

%-Save analysis original max min in seperate files
%--------------------------------------------------------------------------
if (SwE.WB.clusterWise == 1)
  textMaxClusterSize =sprintf('maxClusterSize_WB_%05d', 0);
  if (SwE.WB.stat == 'T')
    textMaxClusterSizeNeg =sprintf('maxClusterSizeNeg_WB_%05d', 0);
  end
end
if (SwE.WB.stat == 'T')
  textMin = sprintf('minScore_WB_%05d', 0);
end
textMax = sprintf('maxScore_WB_%05d', 0);
if spm_matlab_version_chk('7') >=0
  if (SwE.WB.clusterWise == 1)
    save(textMaxClusterSize, 'maxClusterSize','-V6');
    if (SwE.WB.stat == 'T')
      save(textMaxClusterSizeNeg, 'maxClusterSizeNeg','-V6');
    end
  end
  if (SwE.WB.stat == 'T')
    save(textMin,'minScore','-V6');
  end
  save(textMax,'maxScore','-V6');
else
  if (SwE.WB.clusterWise == 1)
    save(textMaxClusterSize, 'maxClusterSize');
    if (SwE.WB.stat == 'T')
      save(textMaxClusterSizeNeg, 'maxClusterSizeNeg');
    end
  end
  if (SwE.WB.stat == 'T')
    save(textMin,'minScore');
  end
  save(textMax,'maxScore');
end

%==========================================================================
%- E N D: Cleanup GUI
%==========================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#
%spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
fprintf('...use the results section for assessment\n\n')
