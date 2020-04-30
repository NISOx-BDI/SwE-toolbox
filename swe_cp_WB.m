function swe_cp_WB(SwE)
% Computes statistic and p-value maps for non-parametric analyses.
% =========================================================================
% For a non-parametric SwE analysis with either NIfTI, GIfTI, CIfTI or '.mat' input, the
% following maps are computed:
%
%   - swe_{unit}_mask: 
%        The mask image for the analysis.
%
%   - swe_{unit}_{T|F}stat_c{c#}: 
%        Voxelwise parametric statistic map (T or F) for contrast {c#}.
%
%   - swe_{unit}_{zT|xF}stat_c{c#}: 
%        Voxelwise parametric equivalent statistic map (Z or Chi Squared) 
%        for contrast {c#}.
%
%   - swe_{unit}_{zT|xF}stat-WB_c{c#}: 
%        Voxelwise non-parametric equivalent statistic map (Z or Chi 
%        Squared) for contrast {c#}.
%
%   - swe_{unit}_{T|F}stat_lp-WB_c{c#}:
%         Log10 map of the voxelwise uncorrected P values for contrast 
%         {c#}.
%
%   - swe_{unit}_{T|F}stat_lpFWE-WB_c{c#}:
%         Log10 map of the voxelwise bootstrap-calculated FWE P values for 
%         contrast {c#}.
%
%   - swe_{unit}_{T|F}stat_lpFDR-WB_c{c#}:
%         Log10 map of the voxelwise bootstrap-calculated FDR P values for
%         contrast {c#}.
%
%   - swe_clustere_{T|F}stat_lpFWE-WB_c{c#}:
%         Log10 map of the clusterwise bootstrap-calculated FWE P values 
%         for contrast {c#}.
%
%   - swe_tfce_c{c#}:
%         TFCE parametric statistic map for contrast {c#}.
%
%   - swe_tfce_lp-WB_c{c#}:
%         Log10 map of the TFCE bootstrap-calculated P values for contrast
%         {c#}.
%    
%   - swe_tfce_lpFWE-WB_c{c#}:
%         Log10 map of the TFCE bootstrap-calculated FWE P values for
%         contrast {c#}.
%
% The field {unit} used above represents the unit in space in which the statistic is calculated. 
% It can be the following strings:
%   - 'vox' for NifTI files,
%   - 'dpx' for GIfTI and CIfTI files,
%   - 'dat' for .mat files.
% Currently (30/08/2018), the only contrasts computed are activation 
% (contrast #1) and deactivation (contrast #2) for the contrast vector the 
% user input during the batch entry.
% =========================================================================
% FORMAT swe_cp_WB(SwE)
% -------------------------------------------------------------------------
% Inputs:
%   - SwE: SwE data structure
% =========================================================================
% Version Info:  $Format:%ci$ $Format:%h$

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

%-Shuffle seed of random number generator
%--------------------------------------------------------------------------
swe_seed

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
file_ext = swe_get_file_extension(SwE.xY.P{1});
isMat    = strcmpi(file_ext,'.mat');
isCifti  = strcmpi(file_ext,'.dtseries.nii') ||  strcmpi(file_ext,'.dscalar.nii');
isOctave = exist('OCTAVE_VERSION','builtin');

if isCifti
  metadata = {'ciftiTemplate', SwE.xY.P{1}};  
  file_data_type = 'dpx';
  dataType = swe_DataType('Cifti');
  dataTypeSpecificInformation = SwE.cifti;
end

if isMat
  file_data_type = 'dat';
  if SwE.WB.clusterWise == 1
    isVolumeMat = isfield(SwE.WB.clusterInfo, 'Vxyz');
    isSurfaceMat = isfield(SwE.WB.clusterInfo, 'Vfaces');
    if isVolumeMat
      dataType = swe_DataType('VolumeMat');
      dataTypeSpecificInformation = [];
    elseif isSurfaceMat
      dataType = swe_DataType('SurfaceMat');
      dataTypeSpecificInformation = importdata(SwE.WB.clusterInfo.Vfaces{1});
      if size(dataTypeSpecificInformation,1) ~=3 && size(dataTypeSpecificInformation,2) ~=3
        error('faces coodinates do not seem correct')
      elseif size(dataTypeSpecificInformation,1) == 3
        dataTypeSpecificInformation = dataTypeSpecificInformation';
      end
    else
      dataType = swe_DataType('Mat');
      dataTypeSpecificInformation = [];
    end
  else
    dataType = swe_DataType('Mat');
    dataTypeSpecificInformation = [];
  end
end

if ~isMat && ~isCifti
  isMeshData = spm_mesh_detect(SwE.xY.VY);
  if isMeshData
      file_ext = '.gii';
      file_data_type = 'dpx';
      dataType = swe_DataType('Gifti');
      g        = SwE.xY.VY(1).private;
      metadata = g.private.metadata;
      name     = {metadata.name};
      if any(ismember(name,'SurfaceID'))
          metadata = metadata(ismember(name,'SurfaceID'));
          metadata = {metadata.name, metadata.value};
      elseif isfield(g,'faces') && ~isempty(g.faces)
          metadata = {'SurfaceID', SwE.xY.VY(1).fname};
      else
          error('SurfaceID not found in GIfTI''s metadata.');
      end
      if isempty(spm_file(metadata{2},'path'))
          metadata{2} = fullfile(spm_file(SwE.xY.VY(1).fname,'path'),metadata{2});
      end
      SwE.xVol.G = metadata{2};
      if (SwE.WB.clusterWise == 1)
          dataTypeSpecificInformation = export(gifti(SwE.xVol.G),'patch');
      end
  else
      dataType = swe_DataType('Nifti');
      dataTypeSpecificInformation = [];
      file_ext = spm_file_ext;
      file_data_type = 'vox';
      metadata = {};
  end
else
  isMeshData = false;
end

try
  giftiAreaFile = SwE.gifti.areaFile;
catch
  giftiAreaFile = '';
end

isVolumeMat = (dataType == swe_DataType('VolumeMat'));
isSurfaceMat = (dataType == swe_DataType('SurfaceMat'));
isNifti = (dataType == swe_DataType('Nifti'));
isGifti  = (dataType == swe_DataType('Gifti'));

%-Check whether we are doing a TFCE analysis
%--------------------------------------------------------------------------
% deactivate for now TFCE if we analyse surface data
TFCE = isfield(SwE.WB, 'TFCE') && ~isMeshData;
if TFCE
    H = SwE.WB.TFCE.H;
    E = SwE.WB.TFCE.E;
    dh = SwE.WB.TFCE.dh;
    C = 18;
end

%-Prevent unnecessary octave warning
%--------------------------------------------------------------------------
if isOctave
   warning ('off', 'histc: empty EDGES specified\n'); 
end

%-Delete files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(SwE.swd,sprintf('swe_v%s_mask%s',file_data_type,file_ext)),'file') == 2
  
  str = {'Current directory contains SwE estimation files:',...
    'pwd = ',SwE.swd,...
    'Existing results will be overwritten!'};
  if spm_input(str,1,'bd','stop|continue',[1,0],1)
    spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
    return
  else
    warning('Overwriting old results\n\t (pwd = %s) ',SwE.swd); %#ok<WNTAG>
  end

end

files = {'^swe_.{3}_mask(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_cov_b\d{2}_b\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_cov_vv(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_cov(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_con_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_beta_b\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_cov_g\d{2}_v\d{2}_v\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_cov_g\d{2}_b\d{2}_b\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_edf_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_beta_\w{1}\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_\w{1,2}stat_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_\w{1,2}stat-WB_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_\w{1,2}stat_lp\w{0,3}_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_\w{1,2}stat_lp\w{0,3}-WB_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_clustere_\w{1,2}stat_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_clustere_\w{1,2}stat_lp\w{0,3}-WB_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_clusternorm\d{0,1}_\w{1,2}stat_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_clusternorm\d{0,1}_\w{1,2}stat_lp\w{0,3}-WB_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_resid_y\d{2,4}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_fit_y\d{2,4}(\.dtseries)?(\.dscalar)?\..{3}$'};

for i = 1:length(files)
  j = spm_select('List',SwE.swd,files{i});
  for k = 1:size(j,1)
    spm_unlink(deblank(j(k,:)));
  end
end

% Tolerance for comparing real numbers
tol = 1e-8;	

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
  if isOctave
    save('SwE','SwE');
  elseif spm_matlab_version_chk('7') >=0
    save('SwE','SwE','-V6');
  else
    save('SwE','SwE');
  end
end

% If clusterWise inference and .mat format, check for the presence of
% spatial information
if isMat && WB.clusterWise == 1
  if isVolumeMat
    XYZ = importdata(SwE.WB.clusterInfo.Vxyz{1});
    if size(XYZ,1) ~=3 && size(XYZ,2) ~=3
      error('voxel coodinates do not seem correct')
    elseif size(XYZ,2) ==3
      XYZ = XYZ';
    end
  elseif ~isSurfaceMat
    error('clusterWise inference cannot be done without spatial information when inputs are in ".mat" format. Please supply faces coordinates (faces or tris) for surface data or voxel coordinates (XYZ_vox) for volumetric data');
  end
end

% small sample correction (for WB)
[corrWB, tmpR2] = swe_resid_corr(SwE, 1, WB.SS, pX);

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
    tmp(i,:) = any(xX.X(iGr_dof==uGr_dof(i),:));
  end
  if nGr_dof==1 | all(sum(tmp, 1)==1) %#ok<OR2>
    break % all is ok, just stop the while
  else
    ind1 = find(sum(tmp, 1)>1,1); % detect the first column in common
    ind2 = find(tmp(:,ind1)==1); % detect the groups to be fused
    for ii = ind2'
      iGr_dof(iGr_dof==uGr_dof(ii)) = ind2(1); % fuse the groups
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
  indSubDesignMatrices = iBeta_dof(ind);
  subjectsInvolved = [];
  for iIndSubDesignMatrices = indSubDesignMatrices
    subjectsInvolved = [subjectsInvolved; iSubj(iGr_dof == iIndSubDesignMatrices)];
  end
  subjectsInvolved = unique(subjectsInvolved);
  indSubjInvolved = nan(length(subjectsInvolved),1);
  % convert into in
  for iSubjInvolved = 1:length(subjectsInvolved)
    indSubjInvolved(iSubjInvolved) = find(uSubj == subjectsInvolved(iSubjInvolved));
  end
  edf = sum(edof_Subj(indSubjInvolved));

  dof_cov = zeros(1,nBeta);
  for i = 1:nBeta
    dof_cov(i) = nSubj_dof(iBeta_dof(i)) - ...
	pB_dof(iBeta_dof(i));    
  end
  
  % This variable should be left empty for Niave estimation.
  dofMat = [];
  xX.erdf_niave = edf;
  SwE.xX = xX;
else
  edf = NaN;
end

%-preprocessing for the modified SwE
if isfield(SwE.type,'modified')
  if dof_type == 1
    error('degrees of freedom type still not implemented for the modified SwE and the WB')
  end
  iVis      = SwE.Vis.iVis;
  iGr       = SwE.Gr.iGr;
  uGr       = unique(iGr); 
  nGr       = length(uGr);
  SwE.Gr.uGr       = uGr;
  SwE.Gr.nGr       = nGr;
    
  % info specific for each group
  uVis_g = cell(1,nGr); % unique visits for each group
  nVis_g = zeros(1,nGr); % number of visits for each group
  uSubj_g = cell(1,nGr); % unique visits for each group
  nSubj_g = zeros(1,nGr); % number of visits for each group
  for g = 1:nGr
    uVis_g{g}  = unique(iVis(iGr==uGr(g))); 
    nVis_g(g)  = length(uVis_g{g});
    iSubj_g = iSubj(iGr==uGr(g)); % Subject number for each subject in group for each visit
    uSubj_g{g} = unique(iSubj_g); % Unique subject numbers of subjects in group
    nSubj_g(g) = length(uSubj_g{g});
    uSubj_g_tmp = uSubj_g{g};
    
    for k = 1:nSubj_g(g)
      
      % The number of visits for subject uSubj_g(k)
      vis_g_subj(k) = sum(iSubj_g==uSubj_g_tmp(k));
      
    end

    max_nVis_g(g) = max(vis_g_subj);
    min_nVis_g(g) = min(vis_g_subj);
        
    clear vis_g_subj
        
  end
  nCov_vis_g  = nVis_g.*(nVis_g+1)/2; % number of covariance elements to be estimated for each group
  nCov_vis    = sum(nCov_vis_g); % total number of covariance elements to be estimated
  
  % Save Vis variables.
  SwE.Vis.uVis_g = uVis_g;
  SwE.Vis.nVis_g = nVis_g;
  SwE.Vis.max_nVis_g = max_nVis_g;
  SwE.Vis.min_nVis_g = min_nVis_g;
  
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
  if dof_type == 1
    edof_Gr = edof_Subj;
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
  DIM      = VY(1).dim;

  YNaNrep = VY(1).dt(2);
    
  fprintf('%-40s: %30s','Output images','...initialising');           %-#
  
  %-Initialise new mask name: current mask & conditions on voxels
  %----------------------------------------------------------------------
  VM    = swe_data_hdr_write(sprintf('swe_%s_mask%s', file_data_type, file_ext), DIM, M,...
                              'swe_cp:resultant analysis mask', metadata, 'uint8');
  
  %-Initialise beta image files
  %----------------------------------------------------------------------

  for i = 1:nBeta
    Vbeta(i) = swe_data_hdr_write(sprintf('swe_%s_beta_b%02d%s',file_data_type,i,file_ext),...
                              DIM, M,...
                              sprintf('swe_cp:beta (%02d) - %s',i,xX.name{i}),...
                              metadata);
  end
  
  %-Initialise original parametric score image, T or F
  %----------------------------------------------------------------------
  if WB.stat=='T'
    eSTAT='z';
  else % F stat
    eSTAT='x';
  end
  
  Vscore = swe_data_hdr_write(sprintf('swe_%s_%cstat_c%02d%s', file_data_type, WB.stat, 1, file_ext), DIM, M,...
			  sprintf('Original parametric %c statistic data.', WB.stat), metadata);
  
  %-Initialise parametric TFCE score image, if TFCE has been selected.
  %---------------------------------------------------------------------- 
  if TFCE
    Vscore_tfce = swe_data_hdr_write(sprintf('swe_tfce_c%02d%s', 1, file_ext), DIM, M,...
				 'Original parametric TFCE statistic data.', metadata); 
    if WB.stat=='T'
      Vscore_tfce_neg = swe_data_hdr_write(sprintf('swe_tfce_c%02d%s', 2, file_ext), DIM, M,...
				       'Original parametric TFCE statistic data for a negative contrast.', metadata); 
    end
  end
  
  %-Initialise original parametric edf image
  %----------------------------------------------------------------------
  
  Vedf = swe_data_hdr_write(sprintf('swe_%s_edf_c%02d%s', file_data_type, 1, file_ext), DIM, M,...
			sprintf('Original parametric %c edf data.', WB.stat), metadata);
          
  %-Initialise parametric P-Value image
  %----------------------------------------------------------------------
  
  VlP = swe_data_hdr_write(sprintf('swe_%s_%cstat_lp_c%02d%s', file_data_type, WB.stat, 1, file_ext), DIM, M,...
		       'Original parametric -log10(P) value data (positive).', metadata);
  
  if WB.stat=='T'
    VlP_Neg = swe_data_hdr_write(sprintf('swe_%s_%cstat_lp_c%02d%s', file_data_type, WB.stat, 2, file_ext), DIM, M,...
			     'Original parametric -log10(P) value data (negative).', metadata);
  end
  
  %-Initialise converted parametric score image
  %----------------------------------------------------------------------
  VcScore = swe_data_hdr_write(sprintf('swe_%s_%c%cstat_c%02d%s', file_data_type, eSTAT, WB.stat, 1, file_ext), DIM, M,...
			   sprintf('Parametric %c statistic data derived from %c-Statistic data.', eSTAT, WB.stat), metadata);    
  if WB.stat=='T'
        VcScore_neg = swe_data_hdr_write(sprintf('swe_%s_%c%cstat_c%02d%s', file_data_type, eSTAT, WB.stat, 2, file_ext), DIM, M,...
			   sprintf('Parametric %c statistic data derived from %c-Statistic data.', eSTAT, WB.stat), metadata);
  end
  
  %-Initialise residual images for the resampling
  %----------------------------------------------------------------------
  
  for i = 1:nScan
    descrip = sprintf('adjusted restricted residuals (%04d)', i);
    VResWB(i) = swe_data_hdr_write(sprintf('swe_%s_resid_y%04d%s', file_data_type, i, file_ext), DIM, M, descrip, metadata);
  end
    
  %-Initialise fitted data images for the resampling
  %----------------------------------------------------------------------
  
  for i = 1:nScan
    descrip = sprintf('restricted fitted data  (%04d)', i);
    VYWB(i) = swe_data_hdr_write(sprintf('swe_%s_fit_y%04d%s',file_data_type,i,file_ext), DIM, M, descrip, metadata);
  end
  
  %-Initialise result images
  %----------------------------------------------------------------------
  VlP_wb_pos = swe_data_hdr_write(sprintf('swe_%s_%cstat_lp-WB_c%02d%s', file_data_type, WB.stat, 1, file_ext), DIM, M,...
                              'Non-parametric voxelwise -log10(P) value data (positive).', metadata);

  VlP_wb_FWE_pos = swe_data_hdr_write(sprintf('swe_%s_%cstat_lpFWE-WB_c%02d%s', file_data_type, WB.stat, 1, file_ext), DIM, M,...
                                  'Non-parametric voxelwise FWE -log10(P) value data (positive).', metadata);
  
  VlP_wb_FDR_pos = swe_data_hdr_write(sprintf('swe_%s_%cstat_lpFDR-WB_c%02d%s', file_data_type, WB.stat, 1, file_ext), DIM, M,...
                                  'Non-parametric voxelwise FDR -log10(P) value data (positive).', metadata);

  if WB.stat=='T'
    VlP_wb_neg = swe_data_hdr_write(sprintf('swe_%s_%cstat_lp-WB_c%02d%s', file_data_type, WB.stat, 2, file_ext), DIM, M,...
				'Non-parametric voxelwise -log10(P) value data (negative).', metadata);
    
    VlP_wb_FWE_neg = swe_data_hdr_write(sprintf('swe_%s_%cstat_lpFWE-WB_c%02d%s', file_data_type, WB.stat, 2, file_ext), DIM, M,...
				    'Non-parametric voxelwise FWE -log10(P) value data (negative).', metadata);
    
    VlP_wb_FDR_neg = swe_data_hdr_write(sprintf('swe_%s_%cstat_lpFDR-WB_c%02d%s', file_data_type, WB.stat, 2, file_ext), DIM, M,...
				    'Non-parametric voxelwise FDR -log10(P) value data (negative).', metadata);

  end
  
  %-Initialise parametric TFCE results images, if TFCE has been selected.
  %---------------------------------------------------------------------- 
  if TFCE
    VlP_tfce_pos = swe_data_hdr_write(sprintf('swe_tfce_lp-WB_c%02d%s', 1, file_ext), DIM, M,...
				  'Non-parametric TFCE -log10(P) value data (positive).', metadata); 
    VlP_tfce_FWE_pos = swe_data_hdr_write(sprintf('swe_tfce_lpFWE-WB_c%02d%s', 1, file_ext), DIM, M,...
				      'Non-parametric TFCE FWE -log10(P) value data (positive).', metadata); 
    if WB.stat=='T'
      VlP_tfce_neg = swe_data_hdr_write(sprintf('swe_tfce_lp-WB_c%02d%s', 2, file_ext), DIM, M,...
				    'Non-parametric TFCE -log10(P) value data (negative).', metadata); 
      VlP_tfce_FWE_neg = swe_data_hdr_write(sprintf('swe_tfce_lpFWE-WB_c%02d%s', 2, file_ext), DIM, M,...
					'Non-parametric TFCE FWE -log10(P) value data (negative).', metadata);
    end
  end
  
  % Converted score for WB.
  VcScore_wb_pos = swe_data_hdr_write(sprintf('swe_%s_%c%cstat-WB_c%02d%s', file_data_type, eSTAT, WB.stat, 1, file_ext), DIM, M,...
                                  sprintf('Non-parametric %c statistic data derived from %c-Statistic data.', eSTAT, WB.stat), metadata);
  
  if WB.clusterWise == 1
      
    % We also need cluster p value maps here.
    V_clustere_pos = swe_data_hdr_write(sprintf('swe_clustere_%cstat_c%02d%s', WB.stat, 1, file_ext), DIM, M,...
                                           sprintf('Cluster extent (positive, CFT %g).',...
                                                   SwE.WB.clusterInfo.primaryThreshold), metadata);
    
    if WB.stat=='T'
      V_clustere_neg = swe_data_hdr_write(sprintf('swe_clustere_%cstat_c%02d%s', WB.stat, 2, file_ext), DIM, M,...
                                          sprintf('Cluster extent (negative, CFT %g).',...
                                                   SwE.WB.clusterInfo.primaryThreshold), metadata);     
    end
  end
  
  fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');    %-#
  %==========================================================================
  % - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
  %==========================================================================
  
  %-Get explicit mask(s)
  %==========================================================================
  mask = true(DIM);
  for i = 1:numel(xM.VM)
    if isCifti
        v = swe_data_read(xM.VM(i)) > 0;
        mask = mask & v(:);
        clear v
    elseif ~(isfield(SwE,'xVol') && isfield(SwE.xVol,'G'))
        %-Assume it fits entirely in memory
        coeff = spm_bsplinc(xM.VM(i), [0 0 0 0 0 0]');
        v = true(DIM);
        [x1,x2] = ndgrid(1:DIM(1),1:DIM(2));
        for x3 = 1:DIM(3)
            M2  = inv(M\xM.VM(i).mat);
            y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*x3+M2(1,4));
            y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*x3+M2(2,4));
            y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*x3+M2(3,4));
            v(:,:,x3) = spm_bsplins(coeff, y1,y2,y3, [0 0 0 0 0 0]') > 0;
        end
        mask = mask & v;
        clear coeff v x1 x2 x3 M2 y1 y2 y3
    else
        if spm_mesh_detect(xM.VM(i))
            v = xM.VM(i).private.cdata() > 0;
        else
            v = spm_mesh_project(gifti(SwE.xVol.G), xM.VM(i)) > 0;
        end
        mask = mask & v(:);
        clear v
    end
  end

  %-Split data into chunks
  %==========================================================================
  chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
  nbchunks  = ceil(prod(DIM) / chunksize);
  chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);
  
  % activated voxels for cluster-wise inference
  if (WB.clusterWise == 1)
    activatedVoxels = false(0);
    maxClusterSize = nan(1, WB.nB + 1);
    if isCifti
      maxClusterSizeInSurfaces = nan(1, WB.nB + 1);
      maxClusterSizeInVolume = nan(1, WB.nB + 1);
    end
    activatedVoxelsNeg = false(0);
    if (WB.stat == 'T')
      maxClusterSizeNeg = nan(1, WB.nB + 1);
      if isCifti
        maxClusterSizeInSurfacesNeg = nan(1, WB.nB + 1);
        maxClusterSizeInVolumeNeg = nan(1, WB.nB + 1);
      end
    end
  end
  maxScore = nan(1, WB.nB + 1);
  if (WB.stat == 'T')
    minScore = nan(1, WB.nB + 1);
  end
  
  %-Cycle over bunches blocks within planes to avoid memory problems
  %==========================================================================
  swe_progress_bar('Init',nbchunks,'Parameter estimation','Chunks');

  for iChunk=1:nbchunks
    chunk = chunks(iChunk):chunks(iChunk+1)-1;

    %-Report progress
    %======================================================================
    if iChunk > 1, fprintf(repmat(sprintf('\b'),1,72)); end                  %-# 
    fprintf('%-40s: %30s', sprintf('Original statistics: Chunk %3d/%-3d',iChunk,nbchunks),...
                              '...processing');
      
    %-Get the data in mask, compute threshold & implicit masks
    %------------------------------------------------------------------
    Y = zeros(nScan, numel(chunk));
    cmask = mask(chunk);
    for iScan=1:nScan
      if ~any(cmask), break, end                 %-Break if empty mask
      
      Y(iScan, cmask) = swe_data_read(VY(iScan), chunk(cmask));%-Read chunk of data
      
      cmask(cmask) = Y(iScan, cmask) > xM.TH(iScan);      %-Threshold (& NaN) mask
      if xM.I && ~YNaNrep && xM.TH(iScan) < 0        %-Use implicit mask
          cmask(cmask) = abs(Y(iScan, cmask)) > eps;
      end
    end
    cmask(cmask) = any(diff(Y(:,cmask),1));  
      
    %-Mask out voxels where data is constant in at least one separable
    % matrix design either in a visit category or within-subject (BG - 27/05/2016)
    %------------------------------------------------------------------
    [cmask, Y, CrS] = swe_mask_seperable(SwE, cmask, Y, iGr_dof);
      
    %==================================================================
    %-Proceed with General Linear Model (if there are voxels)
    %==================================================================
    if CrS
      
      %-General linear model: Ordinary least squares estimation
      %--------------------------------------------------------------        
      beta  = pX*Y;                     %-Parameter estimates
      
      % restricted fitted data
      [resWB, YWB] = swe_fit(SwE, Y, tmpR2, corrWB, beta, SwE.WB.SS);

      if WB.RSwE == 1
        res = swe_fit(SwE, Y, tmpR2, corr, beta, SwE.SS);
      else 
        res = swe_fit(SwE, Y, xX.X, corr, beta, SwE.SS);
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
          cmask(cmask)  = tmp;
          CrS     = sum(cmask);
          Cov_vis = Cov_vis(:,tmp);
        end
        if CrS % Check if there is at least one voxel left
          % compute the visit covariance matrices
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
          if dof_type == 1
            tmpSum = zeros(1,CrS);
          end    
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
              if dof_type == 1
                Cov_beta_g_tmp = weightR(:, iGr_Cov_vis_g==g) * Cov_vis(iGr_Cov_vis_g==g,:);
                if nSizeCon == 1
                  tmpSum = tmpSum + Cov_beta_g_tmp.^2/edof_Gr(g);
                else
                  for iVox = 1:CrS
                    cCovBc_g_vox = zeros(nSizeCon);
                    cCovBc_g_vox(tril(ones(nSizeCon))==1) = Cov_beta_g_tmp(:, iVox);
                    cCovBc_g_vox = cCovBc_g_vox + cCovBc_g_vox' - diag(diag(cCovBc_g_vox));
                    tmpSum(iVox) = tmpSum(iVox) + (trace(cCovBc_g_vox^2) + (trace(cCovBc_g_vox))^2)/edof_Gr(g);
                  end
                end
              end
            end
          end
        end
        cCovBc = weightR * Cov_vis;
        clear Cov_beta_g_tmp
      else % else for "if isfield(SwE.type,'modified')"
        cCovBc = 0;
        if dof_type == 1
          tmpSum = zeros(1,CrS);
        end    
        for i = 1:nSubj
          Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
              (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
          cCovBc = cCovBc + Cov_beta_i_tmp;
          if dof_type == 1
            if nSizeCon == 1 
              tmpSum = tmpSum + Cov_beta_i_tmp.^2/edof_Gr(i);
            else
              for iVox = 1:CrS
                cCovBc_g_vox = zeros(nSizeCon);
                cCovBc_g_vox(tril(ones(nSizeCon))==1) = Cov_beta_i_tmp(:,iVox);
                cCovBc_g_vox = cCovBc_g_vox + cCovBc_g_vox' - diag(diag(cCovBc_g_vox));
                tmpSum(iVox) = tmpSum(iVox) + (trace(cCovBc_g_vox^2) + (trace(cCovBc_g_vox))^2)/edof_Gr(i);
              end
            end
          end
        end
        % These variables are left empty for classic SwE.
        Cov_vis = [];
        dofMat = [];
        clear Cov_beta_i_tmp
      end

      % if dof_type == 1, compute the edf now
      if dof_type == 1
        if nSizeCon == 1
          edf = cCovBc.^2 ./ tmpSum;
        else
          edf = zeros(1,CrS);
          for iVox = 1:CrS
            cCovBc_vox = zeros(nSizeCon);
            cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
            cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
            edf(iVox)=(trace(cCovBc_vox^2) + (trace(cCovBc_vox))^2) / tmpSum(iVox);
          end
        end
      end
      clear tmpSum

      % compute the score
      if (SwE.WB.stat == 'T')
        
        score = (conWB * beta) ./ sqrt(cCovBc);
        
        % hypothesis test, using clusterwise threshold if available.
        if (SwE.WB.clusterWise == 1)
          if dof_type == 1
            hyptest = swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, edf, activatedVoxels, activatedVoxelsNeg);
          else
            hyptest = swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, dofMat, activatedVoxels, activatedVoxelsNeg);
          end
          p = hyptest.positive.p;
          negp = hyptest.negative.p;
          edf = hyptest.positive.edf;
          activatedVoxels = hyptest.positive.activatedVoxels;
          activatedVoxelsNeg = hyptest.negative.activatedVoxels;
          clear CovcCovBc cCovBc
        else
          if dof_type == 1
            hyptest = swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, edf);
          else
            hyptest = swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, dofMat);
          end
          p = hyptest.positive.p;
          negp = hyptest.negative.p;
          edf = hyptest.positive.edf;
        end

        minScore(1) = min(minScore(1), min(hyptest.positive.conScore));

        
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
          if dof_type == 1
            hyptest = swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, edf, activatedVoxels);
          else
            hyptest = swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, dofMat, activatedVoxels);
          end
          p = hyptest.positive.p;
          edf = hyptest.positive.edf;
          activatedVoxels = hyptest.positive.activatedVoxels;
        else
          if dof_type == 1
            hyptest=swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, edf);
          else
            hyptest=swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, dofMat);
          end
          p = hyptest.positive.p;
          edf = hyptest.positive.edf;
        end
      end
      
      maxScore(1) = max(maxScore(1), max(hyptest.positive.conScore));
      
    end % (CrS)

    %-Write output files
    %======================================================================
    c = zeros(numel(chunk),1);

      %-Write mask file
    %----------------------------------------------------------------------
    mask(chunk)  = cmask;
    VM           = swe_data_write(VM, cmask', chunk);
    
    %-Write beta files
    %----------------------------------------------------------------------
    for iBeta = 1:nBeta
      c(cmask) = beta(iBeta,:);
      Vbeta(iBeta) = swe_data_write(Vbeta(iBeta), c, chunk); 
    end

    %-Write WB fitted data images
    %------------------------------------------------------------------
    for iScan = 1:nScan
      c(cmask) = YWB(iScan,:);
      VYWB(iScan) = swe_data_write(VYWB(iScan), c, chunk);
    end
    
    %-Write WB residuals
    %------------------------------------------------------------------
    for iScan = 1:nScan
      c(cmask) = resWB(iScan,:);
      VResWB(iScan) = swe_data_write(VResWB(iScan), c, chunk);
    end

    %-Write parametric score image of the original data
    %------------------------------------------------------------------
    c(cmask) = score;
    Vscore = swe_data_write(Vscore,  c, chunk);
    
    %-Write parametric edf image of the original data
    %------------------------------------------------------------------
    c(cmask) = hyptest.positive.edf;
    Vedf = swe_data_write(Vedf,  c, chunk);
    
    %-Write parametric p-value image
    %------------------------------------------------------------------
    c(cmask) = -log10(hyptest.positive.p);
    VlP = swe_data_write(VlP,  c, chunk);
    
    if WB.stat=='T'
      c(cmask) = -log10(hyptest.negative.p);
      VlP_Neg = swe_data_write(VlP_Neg,  c, chunk);
    end
    
    %-Write converted parametric score image of the original data
    %------------------------------------------------------------------
    c(cmask) = hyptest.positive.conScore;
    VcScore = swe_data_write(VcScore,  c, chunk);

    if WB.stat == 'T'
        VcScore_neg = swe_data_write(VcScore_neg,  -c, chunk);
    end

    %-Report progress
    %======================================================================
    swe_progress_bar('Set',iChunk);
  end % iChunk=1:nbchunks  

  %==========================================================================
  % - P O S T   E S T I M A T I O N   C L E A N U P
  %==========================================================================
  
  S = nnz(mask);
  if S == 0
    error('Please check your data: There are no inmask voxels.');
  end
    
  %-Compute coordinates of voxels within mask
  %--------------------------------------------------------------------------
  [x,y,z]        = ind2sub(DIM,find(mask));
  XYZ            = [x y z]';
    
  if TFCE
    % Create parametric TFCE statistic images.
    if strcmp(WB.stat, 'T')
      
      % Read in T statistics to get negative and positive TFCE scores.
      par_tfce = swe_tfce_transform(swe_data_read(VcScore), H, E, C, dh);
      par_tfce_neg = swe_tfce_transform(-swe_data_read(VcScore), H, E, C, dh);
    else
      
      % Convert F statistics to Z scores.
      scorevol=-swe_invNcdf(10.^(-swe_data_read(VlP)));
      scorevol(isnan(scorevol))=0;
          
      % Convert to TFCE.
      par_tfce = swe_tfce_transform(scorevol, H, E, C, dh);
      
      clear scorevol
    end
    
    % Save parametric TFCE statistic images.
    swe_data_write(Vscore_tfce, par_tfce);
    if strcmp(WB.stat, 'T')
      swe_data_write(Vscore_tfce_neg, par_tfce_neg);
    end
  end
  
  clear beta res Cov_vis c %-Clear to save memory
  
  % compute the max cluster size if needed (so many ways this can be
  % done... Not sure this solution is the best)
  if (SwE.WB.clusterWise == 1)

    LocActivatedVoxels = XYZ(:,activatedVoxels);

    originalClusterStatistics = swe_getClusterStatistics(dataType, LocActivatedVoxels, dataTypeSpecificInformation, giftiAreaFile);
    
    if originalClusterStatistics.nCluster == 0
      warning('no clusters survived the cluster-forming thresholding of the original data for positive effects!')
    end

    maxClusterSize(1) = originalClusterStatistics.maxClusterSize;

    if (SwE.WB.stat == 'T')
      
      LocActivatedVoxelsNeg = XYZ(:,activatedVoxelsNeg);

      originalClusterStatisticsNeg = swe_getClusterStatistics(dataType, LocActivatedVoxelsNeg, dataTypeSpecificInformation, giftiAreaFile);

      if isempty(originalClusterStatisticsNeg.nCluster == 0)
        warning('no clusters survived the cluster-forming thresholding of the original data for negative effects!')
      end

      maxClusterSizeNeg(1) = originalClusterStatisticsNeg.maxClusterSize;

    end
  end

  swe_progress_bar('Clear')

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
  swe_progress_bar('Init',100,str,'');
  
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
  fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...read & mask data')
  
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
  cmask = true(1, nVox);
  %-Use the explicit mask if specified
  if length(SwE.xM.VM) == 1
    cmask(:) = importdata(SwE.xM.VM{1}) > 0;
  end
  %-check if some data need to be masked
  for i = 1:nScan
    if ~any(cmask), break, end                %-Break if empty mask
    cmask(cmask)   = Y(i,cmask) > xM.TH(i);         %-Threshold (& NaN) mask
    if xM.I && ~YNaNrep && xM.TH(i) < 0    %-Use implicit mask
      cmask(cmask) = abs(Y(i,cmask)) > eps;
    end
  end
  %-Mask out voxels where data is constant in at least one separable
  % matrix design either in a visit category or within-subject (BG - 27/05/2016)
  %------------------------------------------------------------------
  [cmask,Y,CrS] = swe_mask_seperable(SwE, cmask, Y, iGr_dof);
  
  if WB.clusterWise == 1
    if isVolumeMat
      XYZ   = XYZ(:,cmask);
    end
  end
  %==================================================================
  %-Proceed with General Linear Model (if there are voxels)
  %==================================================================
  if CrS
    
    %-General linear model: Ordinary least squares estimation
    %--------------------------------------------------------------
    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...estimation');%-#
    
    beta  = pX*Y;                     %-Parameter estimates
    
    [resWB, YWB]=swe_fit(SwE, Y, tmpR2, corrWB, beta, SwE.WB.SS);
    
    if WB.RSwE == 1
      res=swe_fit(SwE, Y, tmpR2, corr, beta, SwE.SS);
    else 
      res=swe_fit(SwE, Y, xX.X, corr, beta, SwE.SS);
    end
    
    clear Y                           %-Clear to save memory
    %-Estimation of the data variance-covariance components (modified SwE)
    %-SwE estimation (classic version)
    %--------------------------------------------------------------
    if dof_type == 1
      tmpSum = zeros(1,CrS);
    end
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
        cmask(cmask)  = tmp;
        CrS     = sum(cmask);
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
          if dof_type == 1
            Cov_beta_g_tmp = weightR(:, iGr_Cov_vis_g==g) * Cov_vis(iGr_Cov_vis_g==g,:);
            if nSizeCon == 1
              tmpSum = tmpSum + Cov_beta_g_tmp.^2/edof_Gr(g);
            else
              for iVox = 1:CrS
                cCovBc_g_vox = zeros(nSizeCon);
                cCovBc_g_vox(tril(ones(nSizeCon))==1) = Cov_beta_g_tmp(:, iVox);
                cCovBc_g_vox = cCovBc_g_vox + cCovBc_g_vox' - diag(diag(cCovBc_g_vox));
                tmpSum(iVox) = tmpSum(iVox) + (trace(cCovBc_g_vox^2) + (trace(cCovBc_g_vox))^2)/edof_Gr(g);
              end
            end
          end
        end
      end
      cCovBc = weightR * Cov_vis;
      clear Cov_beta_g_tmp
    else % classic
      cCovBc = 0;
      for i = 1:nSubj
        Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
	    (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
        cCovBc = cCovBc + Cov_beta_i_tmp;
        if dof_type == 1
          if nSizeCon == 1 
            tmpSum = tmpSum + Cov_beta_i_tmp.^2/edof_Gr(i);
          else
            for iVox = 1:CrS
              cCovBc_g_vox = zeros(nSizeCon);
              cCovBc_g_vox(tril(ones(nSizeCon))==1) = Cov_beta_i_tmp(:,iVox);
              cCovBc_g_vox = cCovBc_g_vox + cCovBc_g_vox' - diag(diag(cCovBc_g_vox));
              tmpSum(iVox) = tmpSum(iVox) + (trace(cCovBc_g_vox^2) + (trace(cCovBc_g_vox))^2)/edof_Gr(i);
            end
          end
        end
      end
      % These variables are left empty for classic SwE.
      Cov_vis = [];
      dofMat = [];
      clear Cov_beta_i_tmp
    end

    % if dof_type == 1, compute the edf now
    if dof_type == 1
      if nSizeCon == 1
        edf = cCovBc.^2 ./ tmpSum;
      else
        edf = zeros(1,CrS);
        for iVox = 1:CrS
          cCovBc_vox = zeros(nSizeCon);
          cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
          cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
          edf(iVox)=(trace(cCovBc_vox^2) + (trace(cCovBc_vox))^2) / tmpSum(iVox);
        end
      end
    end
    clear tmpSum
    % compute the score
    if (SwE.WB.stat == 'T')
      
      score = (conWB * beta) ./ sqrt(cCovBc);
      
      if (SwE.WB.clusterWise == 1)
        if dof_type == 1
          hyptest=swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, edf, activatedVoxels, activatedVoxelsNeg);
        else
          hyptest=swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, dofMat, activatedVoxels, activatedVoxelsNeg);
        end
        p = hyptest.positive.p; 
        negp = hyptest.negative.p;
        edf = hyptest.positive.edf;
        activatedVoxels = hyptest.positive.activatedVoxels;
        activatedVoxelsNeg = hyptest.negative.activatedVoxels;
        clear CovcCovBc cCovBc
      else
        if dof_type == 1
          hyptest=swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, edf);
        else
          hyptest=swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, dofMat);
        end
        p = hyptest.positive.p;
        edf = hyptest.positive.edf;
        clear CovcCovBc cCovBc
      end
         
      minScore(1) = min(hyptest.positive.conScore);
      
      if TFCE
	%%%% TODO: T (MAKE SCORE)
      end
      
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
        if dof_type == 1
          hyptest = swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, edf, activatedVoxels);
        else
          hyptest = swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, dofMat, activatedVoxels);
        end
        p = hyptest.positive.p;
        edf = hyptest.positive.edf;
        activatedVoxels = hyptest.positive.activatedVoxels;
      else
        if dof_type == 1
          hyptest = swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, edf);
        else
          hyptest = swe_hyptest(SwE, score, CrS, cCovBc, Cov_vis, dofMat);
        end
        p = hyptest.positive.p;
        edf = hyptest.positive.edf;
      end
      
      if TFCE
	%%%% TODO: F (MAKE SCORE)
      end
      
    end

    maxScore(1) = max(hyptest.positive.conScore);
    
  end % (CrS)
  M           = [];
  DIM         = [];
  S           = CrS;  
  VM          = sprintf('swe_%s_mask%s', file_data_type, file_ext);
  Vscore      = sprintf('swe_%s_%cstat_c%02d%s', file_data_type, WB.stat, 1, file_ext);
  Vedf        = sprintf('swe_%s_edf_c%02d%s', file_data_type, 1, file_ext);

  mask = cmask;       
  save(sprintf('swe_%s_mask%s', file_data_type, file_ext), 'mask');
  clear mask
  
  edf(:,cmask) = edf;
  save(Vedf, 'edf');
  clear Vedf
  
  tmp = score;
  score = zeros(1, nVox);
  if (SwE.WB.stat == 'T')
    VT(:,cmask) = tmp;
    save(Vscore, 'VT');
    score = tmp;
    clear tmp VT
  else
    VF(:,cmask) = tmp;
    save(Vscore, 'VF');
    score = tmp;
    clear tmp VF
  end
  
  VlP = zeros(1, nVox);
  VlP(:,cmask) = -log10(hyptest.positive.p);
  save(sprintf('swe_%s_%cstat_lp_c%02d%s', file_data_type, WB.stat, 1, file_ext), 'VlP');
  clear VlP
  
  if (SwE.WB.stat == 'T')
    
    VlP_neg = zeros(1, nVox);
    VlP_neg(:,cmask) =  -log10(hyptest.negative.p);
    save(sprintf('swe_%s_%cstat_lp_c%02d%s', file_data_type, WB.stat, 2, file_ext), 'VlP_neg');
    clear VlP_neg
    
    z_map = zeros(1, nVox);
    VZ(:,cmask) =  hyptest.positive.conScore;
    save(sprintf('swe_%s_z%cstat_c%02d%s', file_data_type, WB.stat, 1, file_ext), 'VZ');
    clear VZ
    
  else
      
    x_map = zeros(1, nVox);
    VX(:,cmask) =  hyptest.positive.conScore;
    save(sprintf('swe_%s_x%cstat_c%02d%s', file_data_type, WB.stat, 1, file_ext), 'VX');
    clear VX
    
  end
  
  swe_progress_bar('Clear')

  clear res Cov_vis jj%-Clear to save memory
    
  % compute the max cluster size if needed (so many ways this can be
  % done... Not sure this solution is the best)
  if (SwE.WB.clusterWise == 1)
    
    if dataType == swe_DataType('VolumeMat')
      LocActivatedVoxels = XYZ(:,activatedVoxels);
    else %surface data
      LocActivatedVoxels = false(1, nVox);
      LocActivatedVoxels(cmask) = activatedVoxels;
    end

    originalClusterStatistics = swe_getClusterStatistics(dataType, LocActivatedVoxels, dataTypeSpecificInformation, giftiAreaFile);
    
    if originalClusterStatistics.nCluster == 0
      warning('no clusters survived the cluster-forming thresholding of the original data for positive effects!');
    end

    maxClusterSize(1) = originalClusterStatistics.maxClusterSize;

    if (SwE.WB.stat == 'T')

      if dataType == swe_DataType('VolumeMat')
        LocActivatedVoxelsNeg = XYZ(:,activatedVoxelsNeg);
      else %surface data
        LocActivatedVoxelsNeg = false(nVox,1);
        LocActivatedVoxelsNeg(cmask) = activatedVoxelsNeg;
      end

      originalClusterStatisticsNeg = swe_getClusterStatistics(dataType, LocActivatedVoxelsNeg, dataTypeSpecificInformation, giftiAreaFile);

      if originalClusterStatisticsNeg.nCluster == 0
        warning('no clusters survived the cluster-forming thresholding of the original data for negative effects!')
      end

      maxClusterSizeNeg(1) = originalClusterStatisticsNeg.maxClusterSize;

    end
  end  
end 
%==========================================================================
% - P O S T   E S T I M A T I O N   C L E A N U P
%==========================================================================
if S == 0, spm('alert!','No inmask voxels - empty analysis!'); return; end

%-Save remaining results files and analysis parameters
%==========================================================================

%-place fields in SwE
%--------------------------------------------------------------------------
if WB.clusterWise == 1
  if isfield(SwE.WB, 'clusterInfo') && dataType == swe_DataType('SurfaceMat')
    XYZ = [];
  end
end

if ~isMat
  SwE.xVol.XYZ   = XYZ;               %-InMask XYZ coords (voxels)
end
SwE.xVol.M     = M;                 %-voxels -> mm
SwE.xVol.iM    = inv(M);            %-mm -> voxels
SwE.xVol.DIM   = DIM';               %-image dimensions
SwE.xVol.S     = S;
SwE.xVol.units = {'mm' 'mm' 'mm'};

if ~isMat
  SwE.WB.VYWB       = VYWB;               %-Filehandle - fitted data under H0
  SwE.WB.VResWB     = VResWB;             %-Filehandle - adjusted resticted residuals
  SwE.WB.Vscore     = Vscore;            %-Filehandle - score original image
  SwE.Vbeta      = Vbeta;             %-Filehandle - Beta
end
SwE.WB.weightR    = weightR;
SwE.WB.corrWB     = corrWB;
SwE.WB.corr       = corr;

SwE.WB.tmpR2      = tmpR2;


% cluster-wise specific fields if needed
if (SwE.WB.clusterWise == 1)
  SwE.WB.clusterInfo.LocActivatedVoxels = LocActivatedVoxels;
  SwE.WB.clusterInfo.nCluster = originalClusterStatistics.nCluster;
  SwE.WB.clusterInfo.clusterAssignment = originalClusterStatistics.clusterAssignment;
  SwE.WB.clusterInfo.maxClusterSize = originalClusterStatistics.maxClusterSize;
  SwE.WB.clusterInfo.clusterSize = originalClusterStatistics.clusterSize;
  if isCifti
    SwE.WB.clusterInfo.clusterSizesInSurfaces = originalClusterStatistics.clusterSizesInSurfaces;
    SwE.WB.clusterInfo.clusterSizesInVolume = originalClusterStatistics.clusterSizesInVolume;
  elseif isNifti || isVolumeMat
    SwE.WB.clusterInfo.clusterSizesInVolume = originalClusterStatistics.clusterSize;
  elseif isGifti || isSurfaceMat
    SwE.WB.clusterInfo.clusterSizesInSurfaces = originalClusterStatistics.clusterSize;
  end
  if (SwE.WB.stat == 'T')
    SwE.WB.clusterInfo.LocActivatedVoxelsNeg = LocActivatedVoxelsNeg;
    SwE.WB.clusterInfo.nClusterNeg = originalClusterStatisticsNeg.nCluster;
    SwE.WB.clusterInfo.clusterAssignmentNeg = originalClusterStatisticsNeg.clusterAssignment;
    SwE.WB.clusterInfo.maxClusterSizeNeg = originalClusterStatisticsNeg.maxClusterSize;
    SwE.WB.clusterInfo.clusterSizeNeg = originalClusterStatisticsNeg.clusterSize;
    if isCifti
      SwE.WB.clusterInfo.clusterSizesInSurfacesNeg = originalClusterStatisticsNeg.clusterSizesInSurfaces;
      SwE.WB.clusterInfo.clusterSizesInVolumeNeg = originalClusterStatisticsNeg.clusterSizesInVolume;
    elseif isNifti || isVolumeMat
      SwE.WB.clusterInfo.clusterSizesInVolumeNeg = originalClusterStatisticsNeg.clusterSize;
    elseif isGifti || isSurfaceMat
      SwE.WB.clusterInfo.clusterSizesInSurfacesNeg = originalClusterStatisticsNeg.clusterSize;
    end
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
if isOctave
  save('SwE','SwE');
elseif spm_matlab_version_chk('7') >=0
  save('SwE','SwE','-V6');
else
  save('SwE','SwE');
end

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');             %-#

%==========================================================================
%- Produce bootstraps and maximum stats 
%==========================================================================

% check whether a resampling Matrix exists in SwE. If yes, use it. If not, produce it.
if isfield(SwE.WB, 'resamplingMatrix')
  resamplingMatrix = SwE.WB.resamplingMatrix;
  if any(size(resamplingMatrix) ~= [nScan WB.nB])
    error('The supplied resampling matrix does not have the good dimensions');
  end
else
  % Produce the random value following the Rademacher distribution
  resamplingMatrix = NaN(nScan,WB.nB);
  for iS = 1:nSubj
    % resamplingMatrix(iSubj == uSubj(iS),:) = repmat(binornd(1, 0.5, 1, WB.nB), sum(iSubj == uSubj(iS)), 1);
    resamplingMatrix(iSubj == uSubj(iS),:) = repmat(randi([0 1], 1, WB.nB), sum(iSubj == uSubj(iS)), 1);  % BG (08/11/2016): using randi instead of binornd (which is from the stats toolbox)
  end
  resamplingMatrix(resamplingMatrix == 0) = -1;
end

% load original score
if isMat
  originalScore = hyptest.positive.conScore;
  clear score;
else
  originalScore = swe_data_read(VcScore, 'xyz', XYZ);
  % # blocks
  nbchunks = ceil(S / chunksize);
  chunks = min(cumsum([1 repmat(chunksize, 1, nbchunks)]), S+1);  
end
% variables for results
uncP = ones(1, S); % one because of the original score

str   = sprintf('Parameter estimation\nBootstraping');

swe_progress_bar('Init',100,str,'');

% If we are doing a TFCE analysis we need to record uncorrected P-values
% for TFCE and maxima for FWE.
if TFCE
  tfce_uncP = zeros(DIM(1), DIM(2), DIM(3));
  if SwE.WB.stat == 'T'
    tfce_uncP_neg = zeros(DIM(1), DIM(2), DIM(3));
  end
  
  % We also need to record the TFCE maximas for TFCE FWE (including the
  % original parametric max).
  if TFCE
    maxTFCEScore = nan(1, WB.nB + 1);
    maxTFCEScore(1) = max(par_tfce(:));
    if (WB.stat == 'T')
      maxTFCEScore_neg = nan(1, WB.nB + 1);
      maxTFCEScore_neg(1) = max(par_tfce_neg(:));
    end
  end
  
end

if WB.clusterWise == 1
  clusterSizesInSurfacesUnderH0 = [];
  clusterSizesInVolumeUnderH0 = [];
  if (WB.stat == 'T')
    clusterSizesInSurfacesNegUnderH0 = [];
    clusterSizesInVolumeNegUnderH0 = [];
  end
end

for b = 1:WB.nB
  tic
  str   = sprintf('Parameter estimation\nBootstrap # %i', b);
  swe_progress_bar('Set','xlabel', str)
  
  % activated voxels for cluster-wise inference
  if (SwE.WB.clusterWise == 1)
    activatedVoxels = false(1,S);
    if (SwE.WB.stat == 'T')
      activatedVoxelsNeg = false(1,S);
    end
  end
  if ~isMat
      
    if TFCE
      
      % Instantiate volume for TFCE conversion.
      scorevol = zeros(DIM(1), DIM(2), DIM(3));
      
    end
    
    for iChunk=1:nbchunks
      chunk = chunks(iChunk):chunks(iChunk+1)-1;
      sizeChunk = length(chunk);
      %-Print progress information in command window
      %------------------------------------------------------------------
      if iChunk > 1, fprintf(repmat(sprintf('\b'),1,72)); end                  %-# 
      fprintf('%-40s: %30s', sprintf('Bootstrap # %i: Chunk %3d/%-3d', b, iChunk, nbchunks),...
																	'...processing');
      
      Y_b = swe_data_read(VYWB, 'xyz', XYZ(:,chunk)) + ...
      swe_data_read(VResWB, 'xyz', XYZ(:,chunk)) .* repmat(resamplingMatrix(:,b),1,sizeChunk);
      
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
      if dof_type == 1
        tmpSum = zeros(1,sizeChunk);
      end
      if isfield(SwE.type,'modified')
        Cov_vis=zeros(nCov_vis,sizeChunk);
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
          for iVox = 1:sizeChunk
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
          if dof_type == 1
            Cov_beta_g_tmp = weightR(:, iGr_Cov_vis_g==g) * Cov_vis(iGr_Cov_vis_g==g,:);
            if nSizeCon == 1
              tmpSum = tmpSum + Cov_beta_g_tmp.^2/edof_Gr(g);
            else
              for iVox = 1:sizeChunk
                cCovBc_g_vox = zeros(nSizeCon);
                cCovBc_g_vox(tril(ones(nSizeCon))==1) = Cov_beta_g_tmp(:, iVox);
                cCovBc_g_vox = cCovBc_g_vox + cCovBc_g_vox' - diag(diag(cCovBc_g_vox));
                tmpSum(iVox) = tmpSum(iVox) + (trace(cCovBc_g_vox^2) + (trace(cCovBc_g_vox))^2)/edof_Gr(g);
              end
            end
          end
        end
        cCovBc = weightR * Cov_vis;
      else % classic
        cCovBc = 0;
        for i = 1:nSubj
          Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
	      (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
          cCovBc = cCovBc + Cov_beta_i_tmp;
          if dof_type == 1
            if nSizeCon == 1 
              tmpSum = tmpSum + Cov_beta_i_tmp.^2/edof_Gr(i);
            else
              for iVox = 1:sizeChunk
                cCovBc_g_vox = zeros(nSizeCon);
                cCovBc_g_vox(tril(ones(nSizeCon))==1) = Cov_beta_i_tmp(:,iVox);
                cCovBc_g_vox = cCovBc_g_vox + cCovBc_g_vox' - diag(diag(cCovBc_g_vox));
                tmpSum(iVox) = tmpSum(iVox) + (trace(cCovBc_g_vox^2) + (trace(cCovBc_g_vox))^2)/edof_Gr(i);
              end
            end
          end
        end
        % These variables are left empty for classic SwE.
        Cov_vis = [];
        dofMat = [];
        clear Cov_beta_i_tmp
      end

      % if dof_type == 1, compute the edf now
      if dof_type == 1
        if nSizeCon == 1
          edf = cCovBc.^2 ./ tmpSum;
        else
          edf = zeros(1,sizeChunk);
          for iVox = 1:sizeChunk
            cCovBc_vox = zeros(nSizeCon);
            cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
            cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
            edf(iVox)=(trace(cCovBc_vox^2) + (trace(cCovBc_vox))^2) / tmpSum(iVox);
          end
        end
      end
      clear tmpSum
      % compute the score
      if (SwE.WB.stat == 'T')
	
        score = (conWB * beta) ./ sqrt(cCovBc);
        clear beta
        
      else
	
        cBeta = conWB * beta;
        clear beta
        score = zeros(1, sizeChunk);
        for iVox = 1:sizeChunk
          cCovBc_vox = zeros(nSizeCon);
          cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
          cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
          score(iVox) = cBeta(:,iVox)' / cCovBc_vox * cBeta(:,iVox);
        end
        score = score / rankCon;
        
      end
 
      if dof_type == 1
        hyptest = swe_hyptest(SwE, score, sizeChunk, cCovBc, Cov_vis, edf);
      else
        hyptest = swe_hyptest(SwE, score, sizeChunk, cCovBc, Cov_vis, dofMat);
      end

      if (WB.clusterWise == 1)
        activatedVoxels(chunk) = hyptest.positive.activatedVoxels;
        if (WB.stat == 'T')
          activatedVoxelsNeg(chunk) = hyptest.negative.activatedVoxels;
        end
        clear cCovBc
      end

      uncP(chunk) = uncP(chunk) + (hyptest.positive.conScore > originalScore(chunk) - tol);
      maxScore(b+1) = max(maxScore(b+1), max(hyptest.positive.conScore));
      if (SwE.WB.stat == 'T')
        minScore(b+1) = min(minScore(b+1), min(hyptest.positive.conScore));
      end
      
      % Calculate TFCE uncorrected p image.
      if TFCE

        % Current XYZ indices
        currXYZ = XYZ(1:3, chunk);
	  
        % T test already converted to Z
        if strcmp(WB.stat, 'T')
          scorevol(sub2ind(DIM',currXYZ(1,:),currXYZ(2,:),currXYZ(3,:))) = hyptest.positive.conScore;
        % F test needs to be converted to Z
        else
          % Get score volume from p values
          sv = -swe_invNcdf(hyptest.positive.p);
          % remove NaNs
          sv(isnan(sv))=0;
          % Save as scorevol
          scorevol(sub2ind(DIM',currXYZ(1,:),currXYZ(2,:),currXYZ(3,:))) = sv;
        end
	
      end
      
    end % (iChunk)
    
    if TFCE

      if SwE.WB.stat == 'T'        
	
        % Bootstrapped tfce vol.
        tfce = swe_tfce_transform(scorevol,H,E,C,dh);
        tfce_neg = swe_tfce_transform(-scorevol,H,E,C,dh);    
        
      else
	
        % Bootstrapped tfce vol.
        tfce = swe_tfce_transform(scorevol,H,E,C,dh);
	
      end
      
      % Sum how many voxels are lower than the original parametric tfce.
      tfce_uncP = tfce_uncP + (par_tfce  - tol < tfce);
      if SwE.WB.stat == 'T'
	      tfce_uncP_neg = tfce_uncP_neg + (par_tfce_neg - tol < tfce_neg);
      end
      
      % Record maxima for TFCE FWE p values.
      maxTFCEScore(b+1) = max(tfce(:));
      if SwE.WB.stat == 'T'
	      maxTFCEScore_neg(b+1) = max(tfce_neg(:));
      end
      
      clear tfce tfce_neg
        
    end
    
  else %isMat
    
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
    if dof_type == 1
      tmpSum = zeros(1,S);
    end
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
        if dof_type == 1
          Cov_beta_g_tmp = weightR(:, iGr_Cov_vis_g==g) * Cov_vis(iGr_Cov_vis_g==g,:);
          if nSizeCon == 1
            tmpSum = tmpSum + Cov_beta_g_tmp.^2/edof_Gr(g);
          else
            for iVox = 1:S
              cCovBc_g_vox = zeros(nSizeCon);
              cCovBc_g_vox(tril(ones(nSizeCon))==1) = Cov_beta_g_tmp(:, iVox);
              cCovBc_g_vox = cCovBc_g_vox + cCovBc_g_vox' - diag(diag(cCovBc_g_vox));
              tmpSum(iVox) = tmpSum(iVox) + (trace(cCovBc_g_vox^2) + (trace(cCovBc_g_vox))^2)/edof_Gr(g);
            end
          end
        end
      end
      cCovBc = weightR * Cov_vis;
      clear Cov_beta_g_tmp
    else
      cCovBc = 0;
      for i = 1:nSubj
        Cov_beta_i_tmp = weightR(:,Ind_Cov_vis_classic==i) *...
          (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
        cCovBc = cCovBc + Cov_beta_i_tmp;
        if dof_type == 1
          if nSizeCon == 1 
            tmpSum = tmpSum + Cov_beta_i_tmp.^2/edof_Gr(i);
          else
            for iVox = 1:S
              cCovBc_g_vox = zeros(nSizeCon);
              cCovBc_g_vox(tril(ones(nSizeCon))==1) = Cov_beta_i_tmp(:,iVox);
              cCovBc_g_vox = cCovBc_g_vox + cCovBc_g_vox' - diag(diag(cCovBc_g_vox));
              tmpSum(iVox) = tmpSum(iVox) + (trace(cCovBc_g_vox^2) + (trace(cCovBc_g_vox))^2)/edof_Gr(i);
            end
          end
        end
      end
      
      % These variables are left empty for classic SwE.
      Cov_vis = [];
      dofMat = [];
      clear Cov_beta_i_tmp
    end
 
    % if dof_type == 1, compute the edf now
    if dof_type == 1
      if nSizeCon == 1
        edf = cCovBc.^2 ./ tmpSum;
      else
        edf = zeros(1,S);
        for iVox = 1:S
          cCovBc_vox = zeros(nSizeCon);
          cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
          cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
          edf(iVox)=(trace(cCovBc_vox^2) + (trace(cCovBc_vox))^2) / tmpSum(iVox);
        end
      end
    end
    clear tmpSum

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
    
    if dof_type == 1
      hyptest = swe_hyptest(SwE, score, S, cCovBc, Cov_vis, edf);      
    else
      hyptest = swe_hyptest(SwE, score, S, cCovBc, Cov_vis, dofMat);      
    end

    if (WB.clusterWise == 1)
      activatedVoxels = hyptest.positive.activatedVoxels;
      if (WB.stat == 'T')
          activatedVoxelsNeg = hyptest.negative.activatedVoxels;
      end
      clear cCovBc
    end

    uncP = uncP + (hyptest.positive.conScore > originalScore - tol);
    maxScore(b+1) = max(hyptest.positive.conScore);
    if (SwE.WB.stat == 'T')
      minScore(b+1) = min(hyptest.positive.conScore);
    end     
    
  end
  
  % compute the max cluster size if needed (so many ways this can be
  % done... Not sure this solution is the best)
  if (WB.clusterWise == 1)

    if dataType == swe_DataType('SurfaceMat')
      LocActivatedVoxels = false(1,nVox);
      LocActivatedVoxels(cmask) = activatedVoxels;
    else
      LocActivatedVoxels = XYZ(:,activatedVoxels);
    end

    bootstrapedClusterStatistics = swe_getClusterStatistics(dataType, LocActivatedVoxels, dataTypeSpecificInformation, giftiAreaFile);

    if isfield(bootstrapedClusterStatistics, 'clusterSizesInSurfaces')
      clusterSizesInSurfacesUnderH0 = [clusterSizesInSurfacesUnderH0, bootstrapedClusterStatistics.clusterSizesInSurfaces];
      maxClusterSizeInSurfaces(b+1) = bootstrapedClusterStatistics.maxClusterSizeInSurfaces;
    elseif isfield(bootstrapedClusterStatistics, 'clusterAreas')
      clusterSizesInSurfacesUnderH0 = [clusterSizesInSurfacesUnderH0, bootstrapedClusterStatistics.clusterAreas];
      maxClusterSizeInSurfaces(b+1) = bootstrapedClusterStatistics.maxClusterAreas;
    elseif (dataType == swe_DataType('Gifti') || dataType == swe_DataType('SurfaceMat'))
      clusterSizesInSurfacesUnderH0 = [clusterSizesInSurfacesUnderH0, bootstrapedClusterStatistics.clusterSize];
      maxClusterSizeInSurfaces(b+1) = bootstrapedClusterStatistics.maxClusterSize;
    end
  
    if isfield(bootstrapedClusterStatistics, 'clusterSizesInVolume')
      clusterSizesInVolumeUnderH0 = [clusterSizesInVolumeUnderH0, bootstrapedClusterStatistics.clusterSizesInVolume];
      maxClusterSizeInVolume(b+1) = bootstrapedClusterStatistics.maxClusterSizeInVolume;
    elseif (dataType == swe_DataType('Nifti') || dataType == swe_DataType('VolumeMat'))
      clusterSizesInVolumeUnderH0 = [clusterSizesInVolumeUnderH0, bootstrapedClusterStatistics.clusterSize];
      maxClusterSizeInVolume(b+1) = bootstrapedClusterStatistics.maxClusterSize;
    end

    maxClusterSize(b+1) = bootstrapedClusterStatistics.maxClusterSize;

    if (WB.stat == 'T')

      if dataType == swe_DataType('SurfaceMat')
        LocActivatedVoxelsNeg = false(1,nVox);
        LocActivatedVoxelsNeg(cmask) = activatedVoxelsNeg;
      else
        LocActivatedVoxelsNeg = XYZ(:,activatedVoxelsNeg);
      end

      bootstrapedClusterStatisticsNeg = swe_getClusterStatistics(dataType, LocActivatedVoxelsNeg, dataTypeSpecificInformation, giftiAreaFile);
  
      if isfield(bootstrapedClusterStatisticsNeg, 'maxClusterSizeInVolume')
        maxClusterSizeInVolumeNeg(b+1) = bootstrapedClusterStatisticsNeg.maxClusterSizeInVolume;
      end
  
      if isfield(bootstrapedClusterStatisticsNeg, 'clusterSizesInSurfaces')
        clusterSizesInSurfacesNegUnderH0 = [clusterSizesInSurfacesNegUnderH0, bootstrapedClusterStatisticsNeg.clusterSizesInSurfaces];
        maxClusterSizeInSurfacesNeg(b+1) = bootstrapedClusterStatisticsNeg.maxClusterSizeInSurfaces;
      elseif isfield(bootstrapedClusterStatisticsNeg, 'clusterAreas')
        clusterSizesInSurfacesNegUnderH0 = [clusterSizesInSurfacesNegUnderH0, bootstrapedClusterStatisticsNeg.clusterAreas];
        maxClusterSizeInSurfacesNeg(b+1) = bootstrapedClusterStatisticsNeg.maxClusterArea;
      elseif (dataType == swe_DataType('Gifti') || dataType == swe_DataType('SurfaceMat'))
        clusterSizesInSurfacesNegUnderH0 = [clusterSizesInSurfacesNegUnderH0, bootstrapedClusterStatisticsNeg.clusterSize];
        maxClusterSizeInSurfacesNeg(b+1) = bootstrapedClusterStatisticsNeg.maxClusterSize;
      end
    
      if isfield(bootstrapedClusterStatisticsNeg, 'clusterSizesInVolume')
        clusterSizesInVolumeNegUnderH0 = [clusterSizesInVolumeNegUnderH0, bootstrapedClusterStatisticsNeg.clusterSizesInVolume];
        maxClusterSizeInVolumeNeg(b+1) = bootstrapedClusterStatisticsNeg.maxClusterSizeInVolume;
      elseif (dataType == swe_DataType('Nifti') || dataType == swe_DataType('VolumeMat'))
        clusterSizesInVolumeNegUnderH0 = [clusterSizesInVolumeNegUnderH0, bootstrapedClusterStatisticsNeg.clusterSize];
        maxClusterSizeInVolumeNeg(b+1) = bootstrapedClusterStatisticsNeg.maxClusterSize;
      end
  
      maxClusterSizeNeg(b+1) = bootstrapedClusterStatisticsNeg.maxClusterSize;

    end
  end
  fprintf('%s%30s\n', repmat(sprintf('\b'),1,30), sprintf('..done in %0.4f seconds', toc));
  swe_progress_bar('Set',100 * b / WB.nB);
end

swe_progress_bar('Clear');

%-Save analysis original max min in SwE structure
%--------------------------------------------------------------------------
fprintf('%-40s: %30s','Saving results','...writing');
SwE.WB.maxScore = maxScore;
if (WB.clusterWise == 1)
    SwE.WB.clusterInfo.maxClusterSize = maxClusterSize;
end
if isfield(SwE.WB, 'TFCE')
  SwE.WB.TFCE.maxTFCEScore = maxTFCEScore;
end
if (WB.stat == 'T')
    SwE.WB.minScore = minScore;
    if (WB.clusterWise == 1)
        SwE.WB.clusterInfo.maxClusterSizeNeg = maxClusterSizeNeg;
    end
    if isfield(SwE.WB, 'TFCE')
        SwE.WB.TFCE.maxTFCEScore_neg = maxTFCEScore_neg;
    end
end

%-For cluster-wise analyses, normalise the cluster sizes 
%--------------------------------------------------------------------------
if WB.clusterWise == 1
  scalingFactorNorm = swe_invNcdf(0.75);
  if numel(clusterSizesInSurfacesUnderH0) > 0
    SwE.WB.clusterInfo.clusterSizesInSurfacesUnderH0 = clusterSizesInSurfacesUnderH0;
    lambdaSurfacesUnderH0 = swe_estimateBoxCoxLambda(clusterSizesInSurfacesUnderH0);
    clusterSizesInSurfacesUnderH0 = swe_boxCoxTransform(clusterSizesInSurfacesUnderH0, lambdaSurfacesUnderH0);
    clusterSizesInSurfacesUnderH0_boxCox_mean = mean(clusterSizesInSurfacesUnderH0);
    clusterSizesInSurfacesUnderH0_boxCox_std = std(clusterSizesInSurfacesUnderH0);
    clusterSizesInSurfacesUnderH0_boxCox_median = median(clusterSizesInSurfacesUnderH0);
    clusterSizesInSurfacesUnderH0_boxCox_upperHalfIqr = swe_upperHalfIqr(clusterSizesInSurfacesUnderH0);
    clusterSizesInSurfaces_boxCox = swe_boxCoxTransform(SwE.WB.clusterInfo.clusterSizesInSurfaces, lambdaSurfacesUnderH0);
    
    SwE.WB.clusterInfo.maxClusterSizeInSurfaces = maxClusterSizeInSurfaces;
    SwE.WB.clusterInfo.clusterSizesInSurfacesUnderH0_boxCox_lambda = lambdaSurfacesUnderH0;
    SwE.WB.clusterInfo.clusterSizesInSurfacesUnderH0_boxCox_mean = clusterSizesInSurfacesUnderH0_boxCox_mean;
    SwE.WB.clusterInfo.clusterSizesInSurfacesUnderH0_boxCox_std = clusterSizesInSurfacesUnderH0_boxCox_std;
    SwE.WB.clusterInfo.clusterSizesInSurfacesUnderH0_boxCox_median = clusterSizesInSurfacesUnderH0_boxCox_median;
    SwE.WB.clusterInfo.clusterSizesInSurfacesUnderH0_boxCox_upperHalfIqr = clusterSizesInSurfacesUnderH0_boxCox_upperHalfIqr;
    maxClusterSizeInSurfaces_boxCox = swe_boxCoxTransform(maxClusterSizeInSurfaces, lambdaSurfacesUnderH0);
    
    if (clusterSizesInSurfacesUnderH0_boxCox_upperHalfIqr > 0)
      SwE.WB.clusterInfo.clusterSizesInSurfaces_norm = scalingFactorNorm * ...
        (clusterSizesInSurfaces_boxCox - clusterSizesInSurfacesUnderH0_boxCox_median) ./ clusterSizesInSurfacesUnderH0_boxCox_upperHalfIqr;
      SwE.WB.clusterInfo.maxClusterSizeInSurfaces_norm = scalingFactorNorm * (maxClusterSizeInSurfaces_boxCox - clusterSizesInSurfacesUnderH0_boxCox_median) ./ clusterSizesInSurfacesUnderH0_boxCox_upperHalfIqr;
    else
      SwE.WB.clusterInfo.clusterSizesInSurfaces_norm = nan(1, numel(SwE.WB.clusterInfo.clusterSizesInSurfaces));
      SwE.WB.clusterInfo.maxClusterSizeInSurfaces_norm = nan(1, WB.nB + 1);
    end

  else
    if (isCifti && numel(SwE.cifti.surfaces) > 0) || isMeshData || isSurfaceMat
      warning('no null cluster in surfaces was produced for positive effects!')
      SwE.WB.clusterInfo.clusterSizesInSurfaces_norm = nan(1, numel(SwE.WB.clusterInfo.clusterSizesInSurfaces));
    else
      SwE.WB.clusterInfo.clusterSizesInSurfaces_norm = [];
    end
    SwE.WB.clusterInfo.maxClusterSizeInSurfaces_norm = nan(1, WB.nB + 1);
  end
  
  if numel(clusterSizesInVolumeUnderH0) > 0
    SwE.WB.clusterInfo.clusterSizesInVolumeUnderH0 = clusterSizesInVolumeUnderH0;
    lambdaVolumeUnderH0 = swe_estimateBoxCoxLambda(clusterSizesInVolumeUnderH0);
    clusterSizesInVolumeUnderH0 = swe_boxCoxTransform(clusterSizesInVolumeUnderH0, lambdaVolumeUnderH0);
    clusterSizesInVolumeUnderH0_boxCox_mean = mean(clusterSizesInVolumeUnderH0);
    clusterSizesInVolumeUnderH0_boxCox_std = std(clusterSizesInVolumeUnderH0);
    clusterSizesInVolumeUnderH0_boxCox_median = median(clusterSizesInVolumeUnderH0);
    clusterSizesInVolumeUnderH0_boxCox_upperHalfIqr = swe_upperHalfIqr(clusterSizesInVolumeUnderH0);
    clusterSizesInVolume_boxCox = swe_boxCoxTransform(SwE.WB.clusterInfo.clusterSizesInVolume, lambdaVolumeUnderH0);
    
    SwE.WB.clusterInfo.maxClusterSizeInVolume = maxClusterSizeInVolume;
    SwE.WB.clusterInfo.clusterSizesInVolumeUnderH0_boxCox_lambda = lambdaVolumeUnderH0;
    SwE.WB.clusterInfo.clusterSizesInVolumeUnderH0_boxCox_mean = clusterSizesInVolumeUnderH0_boxCox_mean;
    SwE.WB.clusterInfo.clusterSizesInVolumeUnderH0_boxCox_std = clusterSizesInVolumeUnderH0_boxCox_std;
    SwE.WB.clusterInfo.clusterSizesInVolumeUnderH0_boxCox_median = clusterSizesInVolumeUnderH0_boxCox_median;
    SwE.WB.clusterInfo.clusterSizesInVolumeUnderH0_boxCox_upperHalfIqr = clusterSizesInVolumeUnderH0_boxCox_upperHalfIqr;
    maxClusterSizeInVolume_boxCox = swe_boxCoxTransform(maxClusterSizeInVolume, lambdaVolumeUnderH0);
    
    if (clusterSizesInVolumeUnderH0_boxCox_upperHalfIqr > 0)
      SwE.WB.clusterInfo.clusterSizesInVolume_norm = scalingFactorNorm * ...
        (clusterSizesInVolume_boxCox - clusterSizesInVolumeUnderH0_boxCox_median) ./ clusterSizesInVolumeUnderH0_boxCox_upperHalfIqr;
      SwE.WB.clusterInfo.maxClusterSizeInVolume_norm = scalingFactorNorm * (maxClusterSizeInVolume_boxCox - clusterSizesInVolumeUnderH0_boxCox_median) ./ clusterSizesInVolumeUnderH0_boxCox_upperHalfIqr;
    else
      SwE.WB.clusterInfo.clusterSizesInVolume_norm = nan(1, numel(SwE.WB.clusterInfo.clusterSizesInVolume));
      SwE.WB.clusterInfo.maxClusterSizeInVolume_norm = nan(1, WB.nB + 1);
    end

  else
    if (isCifti && numel(SwE.cifti.volume) > 0) || isNifti || isVolumeMat
      warning('no null cluster in volume was produced for positive effects!')
      SwE.WB.clusterInfo.clusterSizesInVolume_norm = nan(1, numel(SwE.WB.clusterInfo.clusterSizesInVolume));
    else
      SwE.WB.clusterInfo.clusterSizesInVolume_norm = [];
    end
    SwE.WB.clusterInfo.maxClusterSizeInVolume_norm = nan(1, WB.nB + 1);     
  end
  
  SwE.WB.clusterInfo.clusterSize_norm = [SwE.WB.clusterInfo.clusterSizesInSurfaces_norm, SwE.WB.clusterInfo.clusterSizesInVolume_norm];
  SwE.WB.clusterInfo.maxClusterSize_norm = max(SwE.WB.clusterInfo.maxClusterSizeInSurfaces_norm, SwE.WB.clusterInfo.maxClusterSizeInVolume_norm);

  canUseBoxCoxNormalisation = ~any(isnan(SwE.WB.clusterInfo.clusterSize_norm));

  if (WB.stat == 'T')
    if numel(clusterSizesInSurfacesNegUnderH0) > 0
      SwE.WB.clusterInfo.clusterSizesInSurfacesNegUnderH0 = clusterSizesInSurfacesNegUnderH0;
      lambdaSurfacesNegUnderH0 = swe_estimateBoxCoxLambda(clusterSizesInSurfacesNegUnderH0);
      clusterSizesInSurfacesNegUnderH0 = swe_boxCoxTransform(clusterSizesInSurfacesNegUnderH0, lambdaSurfacesNegUnderH0);
      clusterSizesInSurfacesNegUnderH0_boxCox_mean = mean(clusterSizesInSurfacesNegUnderH0);
      clusterSizesInSurfacesNegUnderH0_boxCox_std = std(clusterSizesInSurfacesNegUnderH0);
      clusterSizesInSurfacesNegUnderH0_boxCox_median = median(clusterSizesInSurfacesNegUnderH0);
      clusterSizesInSurfacesNegUnderH0_boxCox_upperHalfIqr = swe_upperHalfIqr(clusterSizesInSurfacesNegUnderH0);
      clusterSizesInSurfacesNeg_boxCox = swe_boxCoxTransform(SwE.WB.clusterInfo.clusterSizesInSurfacesNeg, lambdaSurfacesNegUnderH0);
      
      SwE.WB.clusterInfo.maxClusterSizeInSurfacesNeg = maxClusterSizeInSurfacesNeg;
      SwE.WB.clusterInfo.clusterSizesInSurfacesNegUnderH0_boxCox_lambda = lambdaSurfacesNegUnderH0;
      SwE.WB.clusterInfo.clusterSizesInSurfacesNegUnderH0_boxCox_mean = clusterSizesInSurfacesNegUnderH0_boxCox_mean;
      SwE.WB.clusterInfo.clusterSizesInSurfacesNegUnderH0_boxCox_std = clusterSizesInSurfacesNegUnderH0_boxCox_std;
      SwE.WB.clusterInfo.clusterSizesInSurfacesNegUnderH0_boxCox_median = clusterSizesInSurfacesNegUnderH0_boxCox_median;
      SwE.WB.clusterInfo.clusterSizesInSurfacesNegUnderH0_boxCox_upperHalfIqr = clusterSizesInSurfacesNegUnderH0_boxCox_upperHalfIqr;
      maxClusterSizeInSurfacesNeg_boxCox = swe_boxCoxTransform(maxClusterSizeInSurfacesNeg, lambdaSurfacesNegUnderH0);
      
      if (clusterSizesInSurfacesNegUnderH0_boxCox_upperHalfIqr > 0)
        SwE.WB.clusterInfo.clusterSizesInSurfacesNeg_norm = scalingFactorNorm * ...
          (clusterSizesInSurfacesNeg_boxCox - clusterSizesInSurfacesNegUnderH0_boxCox_median) ./ clusterSizesInSurfacesNegUnderH0_boxCox_upperHalfIqr;
        SwE.WB.clusterInfo.maxClusterSizeInSurfacesNeg_norm = scalingFactorNorm * (maxClusterSizeInSurfacesNeg_boxCox - clusterSizesInSurfacesNegUnderH0_boxCox_median) ./ clusterSizesInSurfacesNegUnderH0_boxCox_upperHalfIqr;
      else
        SwE.WB.clusterInfo.clusterSizesInSurfacesNeg_norm = nan(1, numel(SwE.WB.clusterInfo.clusterSizesInSurfacesNeg));
        SwE.WB.clusterInfo.maxClusterSizeInSurfacesNeg_norm = nan(1, WB.nB + 1);
      end

    else
      if (isCifti && numel(SwE.cifti.surfaces) > 0) || isMeshData || isSurfaceMat
        warning('no null cluster in surfaces was produced for negative effects!')
        SwE.WB.clusterInfo.clusterSizesInSurfacesNeg_norm = nan(1, numel(SwE.WB.clusterInfo.clusterSizesInSurfacesNeg));
      else
        SwE.WB.clusterInfo.clusterSizesInSurfacesNeg_norm = [];
      end
      SwE.WB.clusterInfo.maxClusterSizeInSurfacesNeg_norm = nan(1, WB.nB + 1);
    end

    if numel(clusterSizesInVolumeNegUnderH0) > 0
      SwE.WB.clusterInfo.clusterSizesInVolumeNegUnderH0 = clusterSizesInVolumeNegUnderH0;
      lambdaVolumeNegUnderH0 = swe_estimateBoxCoxLambda(clusterSizesInVolumeNegUnderH0);
      clusterSizesInVolumeNegUnderH0 = swe_boxCoxTransform(clusterSizesInVolumeNegUnderH0, lambdaVolumeNegUnderH0);
      clusterSizesInVolumeNegUnderH0_boxCox_mean = mean(clusterSizesInVolumeNegUnderH0);
      clusterSizesInVolumeNegUnderH0_boxCox_std = std(clusterSizesInVolumeNegUnderH0);
      clusterSizesInVolumeNegUnderH0_boxCox_median = median(clusterSizesInVolumeNegUnderH0);
      clusterSizesInVolumeNegUnderH0_boxCox_upperHalfIqr = swe_upperHalfIqr(clusterSizesInVolumeNegUnderH0);
      clusterSizesInVolumeNeg_boxCox = swe_boxCoxTransform(SwE.WB.clusterInfo.clusterSizesInVolumeNeg, lambdaVolumeNegUnderH0);
      
      SwE.WB.clusterInfo.maxClusterSizeInVolumeNeg = maxClusterSizeInVolumeNeg;
      SwE.WB.clusterInfo.clusterSizesInVolumeNegUnderH0_boxCox_lambda = lambdaVolumeNegUnderH0;
      SwE.WB.clusterInfo.clusterSizesInVolumeNegUnderH0_boxCox_mean = clusterSizesInVolumeNegUnderH0_boxCox_mean;
      SwE.WB.clusterInfo.clusterSizesInVolumeNegUnderH0_boxCox_std = clusterSizesInVolumeNegUnderH0_boxCox_std;
      SwE.WB.clusterInfo.clusterSizesInVolumeNegUnderH0_boxCox_median = clusterSizesInVolumeNegUnderH0_boxCox_median;
      SwE.WB.clusterInfo.clusterSizesInVolumeNegUnderH0_boxCox_upperHalfIqr = clusterSizesInVolumeNegUnderH0_boxCox_upperHalfIqr;
      maxClusterSizeInVolumeNeg_boxCox = swe_boxCoxTransform(maxClusterSizeInVolumeNeg, lambdaVolumeNegUnderH0);
      
      if (clusterSizesInVolumeNegUnderH0_boxCox_upperHalfIqr > 0)
        SwE.WB.clusterInfo.clusterSizesInVolumeNeg_norm = scalingFactorNorm * ...
          (clusterSizesInVolumeNeg_boxCox - clusterSizesInVolumeNegUnderH0_boxCox_median) ./ clusterSizesInVolumeNegUnderH0_boxCox_upperHalfIqr;
        SwE.WB.clusterInfo.maxClusterSizeInVolumeNeg_norm = scalingFactorNorm * (maxClusterSizeInVolumeNeg_boxCox - clusterSizesInVolumeNegUnderH0_boxCox_median) ./ clusterSizesInVolumeNegUnderH0_boxCox_upperHalfIqr;
      else
        SwE.WB.clusterInfo.clusterSizesInVolumeNeg_norm = nan(1, numel(SwE.WB.clusterInfo.clusterSizesInVolumeNeg));
        SwE.WB.clusterInfo.maxClusterSizeInVolumeNeg_norm = nan(1, WB.nB + 1);
      end

    else
      if (isCifti && numel(SwE.cifti.volume) > 0) || isNifti || isVolumeMat
        warning('no null cluster in volume was produced for negative effects!')
        SwE.WB.clusterInfo.clusterSizesInVolumeNeg_norm = nan(1, numel(SwE.WB.clusterInfo.clusterSizesInVolumeNeg));
      else
        SwE.WB.clusterInfo.clusterSizesInVolumeNeg_norm = [];
      end
      SwE.WB.clusterInfo.maxClusterSizeInVolumeNeg_norm = nan(1, WB.nB + 1);
    end

    SwE.WB.clusterInfo.clusterSizeNeg_norm = [SwE.WB.clusterInfo.clusterSizesInSurfacesNeg_norm, SwE.WB.clusterInfo.clusterSizesInVolumeNeg_norm];
    SwE.WB.clusterInfo.maxClusterSizeNeg_norm = max(SwE.WB.clusterInfo.maxClusterSizeInSurfacesNeg_norm, SwE.WB.clusterInfo.maxClusterSizeInVolumeNeg_norm);

    canUseBoxCoxNormalisationNeg = ~any(isnan(SwE.WB.clusterInfo.clusterSizeNeg_norm));

  end
end

%==========================================================================
%- produce results images
%==========================================================================
if isMat
  uncP = uncP / (WB.nB + 1);
  VlP_wb_pos = zeros(1, nVox);
  VlP_wb_pos(:,cmask) = -log10(uncP);
  save(sprintf('swe_%s_%cstat_lp-WB_c%02d%s', file_data_type, WB.stat, 1, file_ext), 'VlP_wb_pos');
  clear VlP_wb_pos
  
  if WB.stat == 'T'
    uncP_neg = 1 + 1/(WB.nB + 1) - uncP;
    VlP_wb_neg = zeros(1, nVox);
    VlP_wb_neg(:,cmask) = -log10(uncP_neg);
    save(sprintf('swe_%s_%cstat_lp-WB_c%02d%s', file_data_type, WB.stat, 2, file_ext), 'VlP_wb_neg');
    clear VlP_wb_neg
    
    VZ_wb = zeros(1, nVox);
    VZ_wb(:,cmask) = swe_invNcdf(1 - uncP);
    save(sprintf('swe_%s_z%cstat-WB_c%02d%s', file_data_type, WB.stat, 1, file_ext), 'VZ_wb');
    clear VZ_wb
    
  else
    
    VX_wb = zeros(1, nVox);
    VX_wb(:,cmask) = spm_invXcdf(1 - uncP,1);
    save(sprintf('swe_%s_x%cstat-WB_c%02d%s', file_data_type, WB.stat, 1, file_ext), 'VX_wb');
    clear VX_wb
    
  end
  
  %
  % - write out lP_FWE+ and lP_FWE- ;
  %
  
  FWERP = ones(1, S); % 1 because the original maxScore is always > original Score
  
  for b = 1:WB.nB
    %-FWER-corrected p is proportion of randomisation greater or
    % equal to statistic.
    %-Use a > b -tol rather than a >= b to avoid comparing
    % two reals for equality.
    FWERP = FWERP + (maxScore(b+1) > originalScore - tol);
  end
  FWERP = FWERP / (WB.nB + 1);
  VlP_wb_FWE_pos = zeros(1, nVox);
  VlP_wb_FWE_pos(:,cmask) = -log10(FWERP);
  save(sprintf('swe_%s_%cstat_lpFWE-WB_c%02d%s',file_data_type,WB.stat,1,file_ext), 'VlP_wb_FWE_pos');
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
    VlP_wb_FWE_neg = zeros(1, nVox);
    VlP_wb_FWE_neg(:,cmask) = -log10(FWERPNeg);
    save(sprintf('swe_%s_%cstat_lpFWE-WB_c%02d%s',file_data_type,WB.stat,2,file_ext), 'VlP_wb_FWE_neg');
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
  VlP_wb_FDR_pos = zeros(1, nVox);
  VlP_wb_FDR_pos(:,cmask) = -log10(fdrP);
  save(sprintf('swe_%s_%cstat_lpFDR-WB_c%02d%s',file_data_type,WB.stat,1,file_ext), 'VlP_wb_FDR_pos');
  clear VlP_wb_FDR_pos fdrP_pos fdrP
  
  if WB.stat =='T'
    try
      fdrP = spm_P_FDR(1 + 1/(WB.nB + 1) - uncP);
    catch
      fdrP = spm_P_FDR(1 + 1/(WB.nB + 1) - uncP,[],'P',[],sort(1 + 1/(WB.nB + 1) - uncP)');
    end
    VlP_wb_FDR_neg = zeros(1, nVox);
    VlP_wb_FDR_neg(:,cmask) = -log10(fdrP);
    save(sprintf('swe_%s_%cstat_lpFDR-WB_c%02d%s',file_data_type,WB.stat,2,file_ext), 'VlP_wb_FDR_neg');
    clear VlP_wb_FDR_neg fdrP_neg fdrP
  end
  
  if WB.clusterWise == 1
    % Not sure what to output. So might be changed later.
    % For now, -log(p_{cluster-wise FWER}) image with nan for non-surviving
    % voxels after the thresholding of the original data
    
    if canUseBoxCoxNormalisation
      
      normClusterFwerP_pos_perCluster = ones(1, SwE.WB.clusterInfo.nCluster); % 1 because the original maxScore is always > original Score
      if (~isempty(SwE.WB.clusterInfo.clusterSize))
        for b = 1:WB.nB
          normClusterFwerP_pos_perCluster = normClusterFwerP_pos_perCluster + (SwE.WB.clusterInfo.maxClusterSize_norm(b+1) > SwE.WB.clusterInfo.clusterSize_norm - tol);
        end
        normClusterFwerP_pos_perCluster = normClusterFwerP_pos_perCluster / (WB.nB + 1);
      end
      
      VlP_wb_normClusterFWE_pos = zeros(1, nVox);
      if isVolumeMat      
        tmp = find(cmask);
        tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxels,2));
        for iC = 1:SwE.WB.clusterInfo.nCluster
          tmp3(SwE.WB.clusterInfo.clusterAssignment == iC) = normClusterFwerP_pos_perCluster(iC);
        end
        VlP_wb_normClusterFWE_pos(tmp(activatedVoxels)) = -log10(tmp3);
      else
        tmp3 = zeros(1, sum(SwE.WB.clusterInfo.LocActivatedVoxels));
        for iC = 1:SwE.WB.clusterInfo.nCluster
          tmp3(SwE.WB.clusterInfo.clusterAssignment == iC) = normClusterFwerP_pos_perCluster(iC);
        end
        VlP_wb_normClusterFWE_pos(SwE.WB.clusterInfo.LocActivatedVoxels) = -log10(tmp3);
      end
      save(sprintf('swe_clusternorm_%cstat_lpFWE-WB_c%02d%s',WB.stat,1,file_ext), 'VlP_wb_normClusterFWE_pos');
    
    else
      
      clusterFwerP_pos_perCluster = ones(1, SwE.WB.clusterInfo.nCluster); % 1 because the original maxScore is always > original Score
      if (~isempty(SwE.WB.clusterInfo.clusterSize))
        for b = 1:WB.nB
          clusterFwerP_pos_perCluster = clusterFwerP_pos_perCluster + (SwE.WB.clusterInfo.maxClusterSize(b+1) > SwE.WB.clusterInfo.clusterSize - tol);
        end
        clusterFwerP_pos_perCluster = clusterFwerP_pos_perCluster / (WB.nB + 1);
      end
      
      VlP_wb_clusterFWE_pos = zeros(1, nVox);
      if isVolumeMat      
        tmp = find(cmask);
        tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxels,2));
        for iC = 1:SwE.WB.clusterInfo.nCluster
          tmp3(SwE.WB.clusterInfo.clusterAssignment == iC) = clusterFwerP_pos_perCluster(iC);
        end
        VlP_wb_clusterFWE_pos(tmp(activatedVoxels)) = -log10(tmp3);
      else
        tmp3 = zeros(1, sum(SwE.WB.clusterInfo.LocActivatedVoxels));
        for iC = 1:SwE.WB.clusterInfo.nCluster
          tmp3(SwE.WB.clusterInfo.clusterAssignment == iC) = clusterFwerP_pos_perCluster(iC);
        end
        VlP_wb_clusterFWE_pos(SwE.WB.clusterInfo.LocActivatedVoxels) = -log10(tmp3);
      end
      save(sprintf('swe_clustere_%cstat_lpFWE-WB_c%02d%s',WB.stat,1,file_ext), 'VlP_wb_clusterFWE_pos');      
    
    end

    if WB.stat =='T'
      
      if canUseBoxCoxNormalisationNeg

        normClusterFwerP_neg_perCluster = ones(1, SwE.WB.clusterInfo.nClusterNeg); % 1 because the original maxScore is always > original Score
        if (~isempty(SwE.WB.clusterInfo.clusterSizeNeg))
          for b = 1:WB.nB
            normClusterFwerP_neg_perCluster = normClusterFwerP_neg_perCluster + (SwE.WB.clusterInfo.maxClusterSizeNeg_norm(b+1) > SwE.WB.clusterInfo.clusterSizeNeg_norm - tol);
          end
          normClusterFwerP_neg_perCluster = normClusterFwerP_neg_perCluster / (WB.nB + 1);
        end
        
        VlP_wb_normClusterFWE_neg = zeros(1, nVox);
        if isVolumeMat
          tmp = find(cmask);
          tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxelsNeg,2));
          for iC = 1:SwE.WB.clusterInfo.nClusterNeg
            tmp3(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = normClusterFwerP_neg_perCluster(iC);
          end
          VlP_wb_normClusterFWE_neg(tmp(activatedVoxelsNeg)) = -log10(tmp3);
        else
          tmp3 = zeros(1, sum(SwE.WB.clusterInfo.LocActivatedVoxelsNeg));
          for iC = 1:SwE.WB.clusterInfo.nClusterNeg
            tmp3(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = normClusterFwerP_neg_perCluster(iC);
          end
          VlP_wb_normClusterFWE_neg(SwE.WB.clusterInfo.LocActivatedVoxelsNeg) = -log10(tmp3);
        end
        save(sprintf('swe_clusternorm_%cstat_lpFWE-WB_c%02d%s',WB.stat,2,file_ext), 'VlP_wb_normClusterFWE_neg');
      
      else
  
        clusterFwerP_neg_perCluster = ones(1, SwE.WB.clusterInfo.nClusterNeg); % 1 because the original maxScore is always > original Score
        if (~isempty(SwE.WB.clusterInfo.clusterSizeNeg))
          for b = 1:WB.nB
            clusterFwerP_neg_perCluster = clusterFwerP_neg_perCluster + (SwE.WB.clusterInfo.maxClusterSizeNeg(b+1) > SwE.WB.clusterInfo.clusterSizeNeg - tol);
          end
          clusterFwerP_neg_perCluster = clusterFwerP_neg_perCluster / (WB.nB + 1);
        end
        
        VlP_wb_clusterFWE_neg = zeros(1, nVox);
        if isVolumeMat
          tmp = find(cmask);
          tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxelsNeg,2));
          for iC = 1:SwE.WB.clusterInfo.nClusterNeg
            tmp3(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = clusterFwerP_neg_perCluster(iC);
          end
          VlP_wb_clusterFWE_neg(tmp(activatedVoxelsNeg)) = -log10(tmp3);
        else
          tmp3 = zeros(1, sum(SwE.WB.clusterInfo.LocActivatedVoxelsNeg));
          for iC = 1:SwE.WB.clusterInfo.nClusterNeg
            tmp3(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = clusterFwerP_neg_perCluster(iC);
          end
          VlP_wb_clusterFWE_neg(SwE.WB.clusterInfo.LocActivatedVoxelsNeg) = -log10(tmp3);
        end
        save(sprintf('swe_clustere_%cstat_lpFWE-WB_c%02d%s',WB.stat,2,file_ext), 'VlP_wb_clusterFWE_neg');
      end
    end
  end
else
  
  Q = cumprod([1,SwE.xVol.DIM(1:2)'])*XYZ - ...
      sum(cumprod(SwE.xVol.DIM(1:2)'));
  %
  % - write out lP+ and lP- images;
  %
  uncP = uncP / (WB.nB + 1);
  tmp = zeros(SwE.xVol.DIM');
  tmp(Q) = -log10(uncP);
  swe_data_write(VlP_wb_pos, tmp);
  
  % If it's F, write out an X map.
  stat = zeros(SwE.xVol.DIM');
  if WB.stat == 'F'
    stat(Q) = spm_invXcdf(1 - uncP,1);
    swe_data_write(VcScore_wb_pos, stat);
  end
  
  % If it's T, write out a Z map.
  if WB.stat == 'T'
      
    % Positive map.
    stat(Q) = swe_invNcdf(1 - uncP);
    swe_data_write(VcScore_wb_pos, stat);

    % T is two tailed so we need a negative map as well.
    tmp(Q) = -log10(1 + 1/(WB.nB + 1) - uncP);
    swe_data_write(VlP_wb_neg, tmp);
  end
  
  %
  % - write out lP_FWE+ and lP_FWE- images;
  %
  
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
  swe_data_write(VlP_wb_FWE_pos, tmp);
  
  % FWE correction for TFCE images.
  if TFCE
      
    % Make new tfce fwe p volume. (Initiate to one to account for
    % original analysis).
    tfcefwevol = ones(DIM(1), DIM(2), DIM(3));
    
    % Calculate FWE p values.
    for b = 1:WB.nB
      tfcefwevol = tfcefwevol + (maxTFCEScore(b+1) > par_tfce - tol);
    end
    tfcefwevol = tfcefwevol / (WB.nB + 1);
      
    % Write out volume.
    swe_data_write(VlP_tfce_FWE_pos, -log10(tfcefwevol));
      
    % Same again for negative contrast, if we are using a T statistic.
    if WB.stat == 'T'
          
      % Make new negative tfce fwe p volume. (Initiate to one to 
      % account for original analysis).
      tfcefwevol_neg = ones(DIM(1), DIM(2), DIM(3));

      % Calculate FWE negative p values.
      for b = 1:WB.nB
	      tfcefwevol_neg = tfcefwevol_neg + (maxTFCEScore_neg(b+1) > par_tfce_neg - tol);
      end
      tfcefwevol_neg = tfcefwevol_neg / (WB.nB + 1);

      % Write out volume.
      swe_data_write(VlP_tfce_FWE_neg, -log10(tfcefwevol_neg));
    end
  end
  
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
    swe_data_write(VlP_wb_FWE_neg, tmp);
  end
  
  %
  % - write out lP_FDR+ and lP_FDR- images;
  %
  try
    tmp(Q) = -log10(spm_P_FDR(uncP));
  catch
    tmp(Q) = -log10(spm_P_FDR(uncP,[],'P',[],sort(uncP)'));
  end
  swe_data_write(VlP_wb_FDR_pos, tmp);
  
  if WB.stat =='T'
    try
      tmp(Q) = -log10(spm_P_FDR(1 + 1/(WB.nB + 1) - uncP));
    catch
      tmp(Q) = -log10(spm_P_FDR(1 + 1/(WB.nB + 1) - uncP,[],'P',[],sort(1 + 1/(WB.nB + 1) - uncP)'));
    end
    swe_data_write(VlP_wb_FDR_neg, tmp);
  end
  
  if WB.clusterWise == 1
    % Not sure what to output. So might be changed later.
    % For now, -log(p_{cluster-wise FWER}) image with nan for non-surviving
    % voxels after the thresholding of the original data
    Q = cumprod([1,SwE.xVol.DIM(1:2)']) * SwE.WB.clusterInfo.LocActivatedVoxels - ...
  	  sum(cumprod(SwE.xVol.DIM(1:2)'));
    tmp= zeros(SwE.xVol.DIM');
    
    tmp4 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxels,2));
    for iC = 1:SwE.WB.clusterInfo.nCluster
      tmp4(SwE.WB.clusterInfo.clusterAssignment == iC) = SwE.WB.clusterInfo.clusterSize(iC);
    end
    tmp(Q) = tmp4;
    swe_data_write(V_clustere_pos, tmp);

    if canUseBoxCoxNormalisation
      
      VlP_wb_normClusterFWE_pos = swe_data_hdr_write(sprintf('swe_clusternorm_%cstat_lpFWE-WB_c%02d%s', WB.stat, 1, file_ext), DIM, M,...
      sprintf('Non-parametric normalised clusterwise FWE -log10(P) value data (positive, CFT %g).',...
      SwE.WB.clusterInfo.primaryThreshold), metadata);

      V_normCluster_pos = swe_data_hdr_write(sprintf('swe_clusternorm_%cstat_c%02d%s', WB.stat, 1, file_ext), DIM, M,...
      sprintf('Box-Cox normalised cluster size (positive, CFT %g).',...
      SwE.WB.clusterInfo.primaryThreshold), metadata);

      normClusterFwerP_pos_perCluster = ones(1, SwE.WB.clusterInfo.nCluster); % 1 because the original maxScore is always > original Score
      if (~isempty(SwE.WB.clusterInfo.clusterSize))
        for b = 1:WB.nB
          normClusterFwerP_pos_perCluster = normClusterFwerP_pos_perCluster + (SwE.WB.clusterInfo.maxClusterSize_norm(b+1) > SwE.WB.clusterInfo.clusterSize_norm - tol);
        end
        normClusterFwerP_pos_perCluster = normClusterFwerP_pos_perCluster / (WB.nB + 1);
      end
      tmp2 = -log10(normClusterFwerP_pos_perCluster);   
      tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxels,2));
      tmp4 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxels,2));
      for iC = 1:SwE.WB.clusterInfo.nCluster
        tmp3(SwE.WB.clusterInfo.clusterAssignment == iC) = tmp2(iC);
        tmp4(SwE.WB.clusterInfo.clusterAssignment == iC) = SwE.WB.clusterInfo.clusterSize_norm(iC);
      end
      tmp(Q) = tmp3;
      swe_data_write(VlP_wb_normClusterFWE_pos, tmp);
      tmp(Q) = tmp4;
      swe_data_write(V_normCluster_pos, tmp);
      
    else

      VlP_wb_clusterFWE_pos = swe_data_hdr_write(sprintf('swe_clustere_%cstat_lpFWE-WB_c%02d%s', WB.stat, 1, file_ext), DIM, M,...
      sprintf('Non-parametric clusterwise FWE -log10(P) value data (positive, CFT %g).',...
      SwE.WB.clusterInfo.primaryThreshold), metadata);

      clusterFwerP_pos_perCluster = ones(1, SwE.WB.clusterInfo.nCluster); % 1 because the original maxScore is always > original Score
      if (~isempty(SwE.WB.clusterInfo.clusterSize))
        for b = 1:WB.nB
          clusterFwerP_pos_perCluster = clusterFwerP_pos_perCluster + (SwE.WB.clusterInfo.maxClusterSize(b+1) > SwE.WB.clusterInfo.clusterSize - tol);
        end
        clusterFwerP_pos_perCluster = clusterFwerP_pos_perCluster / (WB.nB + 1);
      end
      tmp2 = -log10(clusterFwerP_pos_perCluster);   
      tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxels,2));
      for iC = 1:SwE.WB.clusterInfo.nCluster
        tmp3(SwE.WB.clusterInfo.clusterAssignment == iC) = tmp2(iC);
      end
      tmp(Q) = tmp3;
      swe_data_write(VlP_wb_clusterFWE_pos, tmp);

    end
    
    if WB.stat =='T'
      Q = cumprod([1,SwE.xVol.DIM(1:2)']) * SwE.WB.clusterInfo.LocActivatedVoxelsNeg - ...
        sum(cumprod(SwE.xVol.DIM(1:2)'));
      tmp= zeros(SwE.xVol.DIM');
      
      tmp4 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxelsNeg, 2));
      for iC = 1:SwE.WB.clusterInfo.nClusterNeg
        tmp4(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = SwE.WB.clusterInfo.clusterSizeNeg(iC);
      end
      tmp(Q) = tmp4;
      swe_data_write(V_clustere_neg, tmp);

      if canUseBoxCoxNormalisationNeg

        VlP_wb_normClusterFWE_neg = swe_data_hdr_write(sprintf('swe_clusternorm_%cstat_lpFWE-WB_c%02d%s', WB.stat, 2, file_ext), DIM, M,...
        sprintf('Non-parametric normalised clusterwise FWE -log10(P) value data (negative, CFT %g).',...
        SwE.WB.clusterInfo.primaryThreshold), metadata);

        V_normCluster_neg = swe_data_hdr_write(sprintf('swe_clusternorm_%cstat_c%02d%s', WB.stat, 2, file_ext), DIM, M,...
        sprintf('Box-Cox normalised cluster size (negative, CFT %g).',...
        SwE.WB.clusterInfo.primaryThreshold), metadata);

        normClusterFwerP_neg_perCluster = ones(1, SwE.WB.clusterInfo.nClusterNeg); % 1 because the original maxScore is always > original Score
        if (~isempty(SwE.WB.clusterInfo.clusterSizeNeg))
          for b = 1:WB.nB
            normClusterFwerP_neg_perCluster = normClusterFwerP_neg_perCluster + (SwE.WB.clusterInfo.maxClusterSizeNeg_norm(b+1) > SwE.WB.clusterInfo.clusterSizeNeg_norm - tol);
          end
          normClusterFwerP_neg_perCluster = normClusterFwerP_neg_perCluster / (WB.nB + 1);
        end
        tmp2 = -log10(normClusterFwerP_neg_perCluster);
        tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxelsNeg,2));
        tmp4 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxelsNeg,2));
        for iC = 1:SwE.WB.clusterInfo.nClusterNeg
          tmp3(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = tmp2(iC);
          tmp4(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = SwE.WB.clusterInfo.clusterSizeNeg_norm(iC);
        end
        tmp(Q) = tmp3;
        swe_data_write(VlP_wb_normClusterFWE_neg, tmp);
        tmp(Q) = tmp4;
        swe_data_write(V_normCluster_neg, tmp);

      else

        VlP_wb_clusterFWE_neg = swe_data_hdr_write(sprintf('swe_clustere_%cstat_lpFWE-WB_c%02d%s', WB.stat, 2, file_ext), DIM, M,...
        sprintf('Non-parametric clusterwise FWE -log10(P) value data (negative, CFT %g).',...
        SwE.WB.clusterInfo.primaryThreshold), metadata);

        clusterFwerP_neg_perCluster = ones(1, SwE.WB.clusterInfo.nClusterNeg); % 1 because the original maxScore is always > original Score
        if (~isempty(SwE.WB.clusterInfo.clusterSizeNeg))
          for b = 1:WB.nB
            clusterFwerP_neg_perCluster = clusterFwerP_neg_perCluster + (SwE.WB.clusterInfo.maxClusterSizeNeg(b+1) > SwE.WB.clusterInfo.clusterSizeNeg - tol);
          end
          clusterFwerP_neg_perCluster = clusterFwerP_neg_perCluster / (WB.nB + 1);
        end
        tmp2 = -log10(clusterFwerP_neg_perCluster);
        tmp3 = zeros(1, size(SwE.WB.clusterInfo.LocActivatedVoxelsNeg,2));
        for iC = 1:SwE.WB.clusterInfo.nClusterNeg
          tmp3(SwE.WB.clusterInfo.clusterAssignmentNeg == iC) = tmp2(iC);
        end
        tmp(Q) = tmp3;
        swe_data_write(VlP_wb_clusterFWE_neg, tmp);

      end
      
    end
  end
  
  if TFCE
    tfce_luncP = -log10((tfce_uncP+1)./(WB.nB+1));
    swe_data_write(VlP_tfce_pos, tfce_luncP);
    if WB.stat == 'T'
      tfce_luncP_neg = -log10((tfce_uncP_neg+1)./(WB.nB+1));
      swe_data_write(VlP_tfce_neg, tfce_luncP_neg);
    end
  end
  
end

% save the version number of the toolbox
SwE.ver = swe('ver');

SwE.xY.dataType = dataType;

%==========================================================================
%- E N D: Cleanup GUI
%==========================================================================

if ~isMat
  % Remove residual and Y images now we are done with them:
  files = {'^swe_.{3}_resid_y\d{2,4}(\.dtseries)?(\.dscalar)?\..{3}$','^swe_.{3}_fit_y\d{2,4}(\.dtseries)?(\.dscalar)?\..{3}$'};
  for i = 1:numel(files)
    j = cellstr(spm_select('FPList',SwE.swd,files{i}));
    for k = 1:numel(j)
      spm_unlink(j{k});
    end
  end
end

%Save SwE.
if isOctave
  save('SwE','SwE');
elseif spm_matlab_version_chk('7') >=0
  save('SwE','SwE','-V6');
else
  save('SwE','SwE');
end

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#
%spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
fprintf('...use the saved images for assessment\n\n')

end

%-Mask out voxels where data is constant in at least one separable
% matrix design either in a visit category or within-subject (BG - 27/05/2016)
function [cmask,Y,CrS]=swe_mask_seperable(SwE, cmask, Y, iGr_dof)
    
% Setup
nGr_dof = length(unique(iGr_dof));
if isfield(SwE.type,'modified')
  nGr = SwE.Gr.nGr;
  iGr = SwE.Gr.iGr;
  uGr = SwE.Gr.uGr;
  iVis = SwE.Vis.iVis;
  iSubj = SwE.Subj.iSubj;
  nVis_g = SwE.Vis.nVis_g;
  uVis_g = SwE.Vis.uVis_g;
end

% first look data for each separable matrix design
for g = 1:nGr_dof
  % do not look for cases where the separable matrix design is only one row (BG - 05/08/2016)
  if sum(iGr_dof'==g) > 1
    % mask constant data within separable matrix design g (added by BG on 29/08/16)
    cmask(cmask) = any(abs(diff(Y(iGr_dof'==g,cmask),1)) > eps, 1);
    if isfield(SwE.type,'modified') % added by BG on 29/08/16
      % then look data for each "homogeneous" group check if the data is contant over subject for each visit category
      for g2 = 1:nGr
        for k = 1:nVis_g(g2)
          if sum(iGr_dof'==g & iGr == uGr(g2) & iVis == uVis_g{g2}(k)) > 1 % do not look for cases when the data is only one row (BG - 05/08/2016)
            cmask(cmask) = any(abs(diff(Y(iGr_dof'==g & iGr == uGr(g2) & iVis == uVis_g{g2}(k) ,cmask),1)) > eps, 1);
            for kk = k:nVis_g(g2)
              if k ~= kk
                % extract the list of subject with both visit k and kk
                subjList = intersect(iSubj(iGr_dof'==g & iGr == uGr(g2) & iVis == uVis_g{g2}(k)), iSubj(iGr_dof'==g & iGr == uGr(g2) & iVis == uVis_g{g2}(kk)));
                % look if some difference are observed within subject
                if ~isempty(subjList)
                  diffVis = cmask(cmask) == 0;
                  for i = 1:length(subjList)
                    diffVis = diffVis | (abs(Y(iSubj == subjList(i) & iVis == uVis_g{g2}(k), cmask) - Y(iSubj == subjList(i) & iVis == uVis_g{g2}(kk), cmask)) > eps);
                  end
                  cmask(cmask) = diffVis;
                end
              end
            end
          end
        end
      end
    end
  end
end
      
Y      = Y(:,cmask);                          %-Data within mask
CrS    = sum(cmask);                          %-# current voxels

end

% This function performs a hypothesis test using the threshold given as the
% primary threshold in the SwE cluster info. If this is not available it
% returns unthresholded values only.
function hyptest =swe_hyptest(SwE, score, matSize, cCovBc, Cov_vis, dofMat, varargin)

% setup
p = zeros(1, matSize);
nSizeCon = size(SwE.WB.con,1);
rankCon = rank(SwE.WB.con);

if isfield(SwE.type,'modified')
  dof_type = SwE.type.modified.dof_mo;
  nGr = length(unique(SwE.Gr.iGr));
  
  if nSizeCon == 1
    Wg_2 = SwE.WB.Wg{1};
    Wg_3 = SwE.WB.Wg{1};
  else
    Wg_2 = SwE.WB.Wg{2};
    Wg_3 = SwE.WB.Wg{3};
  end
  
else
  dof_type = SwE.type.classic.dof_cl;
  nGr = SwE.Subj.nSubj;        
end

% Convert P values.
switch dof_type
    
  case 0
      
    edf = SwE.xX.erdf_niave;
    
  case 1
    % for dof_type = 1, saved in dofMat. Should be changed later
    edf = dofMat;
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
  % We calculate both p and negative p as both will experience overflow
  % around 1 which may be undesirable if the user wishes to view both
  % results.
  p  = spm_Tcdf(-score, edf);
  negp = spm_Tcdf(score, edf);
  
  if SwE.WB.clusterWise~=0
    if nargin <=6
      % We may wish to just record the activated voxels. 
      activatedVoxels = p < (SwE.WB.clusterInfo.primaryThreshold);
      activatedVoxelsNeg = negp < (SwE.WB.clusterInfo.primaryThreshold);
    else
      % Or we may wish to add the activatedVoxels to a pre-existing list.
      activatedVoxels = [varargin{1}, p < (SwE.WB.clusterInfo.primaryThreshold)];
      activatedVoxelsNeg = [varargin{2}, negp < (SwE.WB.clusterInfo.primaryThreshold)];
    end
  end
  
else
  scoreTmp = (edf-rankCon+1) ./ edf .* score;
  scoreTmp(scoreTmp < 0 ) = 0;
  if dof_type == 0
    p(scoreTmp>0) = betainc((edf-rankCon+1)./(edf-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf-rankCon+1)/2, rankCon/2);
  else
    p(scoreTmp>0) = betainc((edf(scoreTmp>0)-rankCon+1)./(edf(scoreTmp>0)-rankCon+1+rankCon*scoreTmp(scoreTmp>0)),(edf(scoreTmp>0)-rankCon+1)/2, rankCon/2);
    p(scoreTmp <= 0) = 1;
  end

  if SwE.WB.clusterWise~=0
    if nargin<=6
      activatedVoxels = p < (SwE.WB.clusterInfo.primaryThreshold);
    else
      activatedVoxels = [varargin{1}, p < (SwE.WB.clusterInfo.primaryThreshold)];
    end
  end
  
  negp = 1-p;
  activatedVoxelsNeg = NaN;
  
end

% Calculate converted score in a way which minimizes under/overflow
if (SwE.WB.stat == 'T')
    
    conScore = zeros(size(score));
    
    switch dof_type
        % In case 0 edf is the same for everything
        case 0
            
            if any(score>0)
                conScore(score>0) = -spm_invNcdf(spm_Tcdf(-score(score>0),edf));
            end
            if any(score<=0)
                conScore(score<=0) = spm_invNcdf(spm_Tcdf(score(score<=0),edf));
            end
            
        % In all other cases edf varies over voxels
        otherwise
            
            if any(score>0)
                conScore(score>0) = -spm_invNcdf(spm_Tcdf(-score(score>0),edf(score>0)));
            end
            if any(score<=0)
                conScore(score<=0) = spm_invNcdf(spm_Tcdf(score(score<=0),edf(score<=0)));
            end

    end
    
else
    conScore = swe_invNcdf(0.5 * p).^2;
end
    
% Save results
hyptest = struct('positive', struct('p', p,...
                                    'edf', edf,...
                                    'conScore', conScore));

% Save activated voxels if doing clusterwise.
if SwE.WB.clusterWise~=0
    hyptest.positive.activatedVoxels = activatedVoxels;
end
                                   
% If it's a T contrast save negative results as well.
if (SwE.WB.stat == 'T')
    hyptest.negative = struct('p', negp,...
                              'edf', edf,...
                              'conScore', -conScore);

    % Save activated voxels if doing clusterwise.
    if SwE.WB.clusterWise~=0
        hyptest.negative.activatedVoxels = activatedVoxelsNeg;
    end

else
    
    % For an F, negative p-values are still useful
    hyptest.negative = struct('p', negp);
    
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
    if restric == 1
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

function value = swe_upperHalfIqr(X)
  value = diff( quantile(X, [.5 .75]) );
end