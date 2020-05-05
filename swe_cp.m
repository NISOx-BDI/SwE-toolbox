function swe_cp(SwE)
% This function computes covariance and beta maps for parametric analyses.
% =========================================================================
% For a parametric SwE analysis with nifti input, this function computes 
% the following maps:
%
%   - swe_vox_mask: The mask image for the analysis.
%   - swe_vox_con_c{c#}: The contrast map for contrast {c#}
%   - swe_vox_cov_b{b1#}_b{b2#}: The covariance map between betas {b1#}
%     and {b2#}.
%   - swe_vox_cov_g{g#}_b{b1#}_b{b2#}: The covariance map between betas
%     {b1#} and {b2#} for group {g#}.
%   - swe_vox_cov_g{g#}_v{v1#}_v{v2#}: The covariance map between betas
%     {v1#} and {v2#} for group {g#}.
%
% For a parametric SwE analysis with GIfTI or CIfTI inputs, this function computes 
% the following maps:
%
%   - swe_dpx_mask: The mask image for the analysis.
%   - swe_dpx_con_c{c#}: The contrast map for contrast {c#}
%   - swe_dpx_cov_b{b1#}_b{b2#}: The covariance map between betas {b1#}
%     and {b2#}.
%   - swe_dpx_cov_g{g#}_b{b1#}_b{b2#}: The covariance map between betas
%     {b1#} and {b2#} for group {g#}.
%   - swe_dpx_cov_g{g#}_v{v1#}_v{v2#}: The covariance map between betas
%     {v1#} and {v2#} for group {g#}.
%
% For a parametric SwE analysis with '.mat' input, this function computes
% the following analagous maps:
%
%   - swe_dat_mask: The mask image for the analysis.
%   - swe_dat_beta_b: The beta map.
%   - swe_dat_con_c: The contrast map for each contrast.
%   - swe_dat_cov_bb: The between-betas covariance map.
%   - swe_dat_cov_g_bb: The groupwise between-betas covariance maps.
%   - swe_dat_cov_g_vv: The visitwise between-betas covariance maps.
%
% For non-parametric SwE analyses, the function `swe_cp_WB` is called
% instead as these maps must be computed differently. See the header of
% `swe_cp_WB` for more information.
% =========================================================================
% FORMAT swe_cp(SwE)
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

%-Get SwE.mat[s] if necessary
%--------------------------------------------------------------------------
if nargin == 0
    P     = cellstr(spm_select(Inf,'^SwE\.mat$','Select SwE.mat[s]'));
    for i = 1:length(P)
        swd     = fileparts(P{i});
        load(fullfile(swd,'SwE.mat'));
        SwE.swd = swd;
        
        % detect if this is a WB analysis or a "standard analysis"
        if isfield(SwE, 'WB')
          swe_cp_WB(SwE);
        else
          swe_cp(SwE);
        end
    end
    return
end
% If this is a WB analysis we need to use swe_cp_WB.
if isfield(SwE, 'WB')
     swe_cp_WB(SwE);
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

%-Check if we have data in a.mat format and set some variables accordingly
%--------------------------------------------------------------------------
file_ext = swe_get_file_extension(SwE.xY.P{1});
isMat    = strcmpi(file_ext,'.mat');
isCifti  = strcmpi(file_ext,'.dtseries.nii') ||  strcmpi(file_ext,'.dscalar.nii');
isOctave = exist('OCTAVE_VERSION','builtin');

if isCifti
    metadata = {'ciftiTemplate', SwE.xY.P{1}};
    file_data_type = 'dpx';
end

if isMat
    file_data_type = 'dat';
end

if ~isMat && ~isCifti
    isMeshData = spm_mesh_detect(SwE.xY.VY);
    if isMeshData
        file_ext = '.gii';
        file_data_type = 'dpx';
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
    else
        file_ext = spm_file_ext;
        file_data_type = 'vox';
        metadata = {};
    end
end

%-Delete files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(SwE.swd,sprintf('swe_%s_mask%s',file_data_type,file_ext)),'file') == 2
 
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
        '^swe_clustere_\w{1,2}stat_lp\w{0,3}_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_clustere_\w{1,2}stat_lp\w{0,3}-WB_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_clusternorm_\w{1,2}stat_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_clusternorm_\w{1,2}stat_lp\w{0,3}-WB_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_clusternorm2_\w{1,2}stat_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_clusternorm2_\w{1,2}stat_lp\w{0,3}-WB_c\d{2}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_resid_y\d{2,4}(\.dtseries)?(\.dscalar)?\..{3}$',...
        '^swe_.{3}_fit_y\d{2,4}(\.dtseries)?(\.dscalar)?\..{3}$'};
 
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
Hat           = xX.X*(pX); % Hat matrix
iSubj         = SwE.Subj.iSubj;
uSubj         = unique(iSubj);
nSubj         = length(uSubj);

%-residual correction
%
switch SwE.SS
    case 0
        corr = ones(nScan,1);
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
    if nGr_dof==1 | all(sum(tmp,1)==1) %#ok<OR2>
        break % all is ok, just stop the while
    else
        ind1 = find(sum(tmp,1)>1,1); % detect the first column in common
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
    % weights for the vectorised SwE (to be checked)
    weight=NaN(nCov_beta,nCov_vis);
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
    it = 0;    
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
            [~,n,e]     = fileparts(VY(i).fname);
            VY(i).fname = [n,e];
        end
    end

    M        = VY(1).mat;
    DIM      = VY(1).dim;

    % check how the data image treat 0 (as NaN or not)
    YNaNrep = VY(1).dt(2);

    %-Maximum number of residual images for smoothness estimation
    %--------------------------------------------------------------------------
    % MAXRES   = Inf; (commented by BG on 08/11/2016)
    % nSres    = nScan; (commented by BG on 08/11/2016)

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

    %-Initialise Cov_beta image files
    %----------------------------------------------------------------------

    it=0;
    for i=1:nBeta
        for ii=i:nBeta
            it=it+1;
            Vcov_beta(it) = swe_data_hdr_write(sprintf('swe_%s_cov_b%02d_b%02d%s',file_data_type,i,ii,file_ext),...
                                           DIM, M, sprintf('cov_beta_%02d_%02d hats - %s/%s',...
                                                i,ii,xX.name{i},xX.name{ii}),...
                                           metadata);
        end
    end


    %-Initialise Cov_beta_g image files if needed
    %----------------------------------------------------------------------
    if dof_type == 1
        if isfield(SwE.type,'classic')
            nGr = nSubj;
            uGr = uSubj;
        end

        it=0;
        for g=1:nGr
            for ii=1:nBeta
                for iii=ii:nBeta
                    it=it+1;
                    Vcov_beta_g(it) = swe_data_hdr_write([sprintf('swe_%s_cov_g%02d_b%02d_b%02d',file_data_type,g,ii,iii) file_ext],...
                        DIM, M, sprintf('cov_beta_g_%02d_%02d_%02d hats - group %s - %s/%s',...
                            g,ii,iii,num2str(uGr(g)),xX.name{ii},xX.name{iii}), metadata);
                end
            end
        end
    end

    %-Initialise cov_vis image files
    %----------------------------------------------------------------------
    if isfield(SwE.type,'modified')
        it=0;
        for g =1:nGr
            for ii=1:nVis_g(g)
                for iii=ii:nVis_g(g)
                    it=it+1;
                    Vcov_vis(it) = swe_data_hdr_write([sprintf('swe_%s_cov_g%02d_v%02d_v%02d',file_data_type,g,ii,iii) file_ext],...
                                                  DIM, M, sprintf('cov_vis_%02d_%02d_%02d hats - group %s - visits %s/%s',...
                                                       g,ii,iii,num2str(uGr(g)),num2str(uVis_g{g}(ii)),num2str(uVis_g{g}(iii))),...
                                                       metadata);
                end
            end
        end
    end
    %-Initialise standardised residual images
    %----------------------------------------------------------------------
    % for i = 1:nSres
    %     VResI(i) = swe_create_vol(sprintf('swe_%s_resid_y%02d.img', i),...
    %                               DIM, M, sprintf('spm_spm:ResI (%02d)', i),...
    %                               isMeshData);
    % end
    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');    %-# 
    %  
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
          C = spm_bsplinc(xM.VM(i), [0 0 0 0 0 0]');
          v = true(DIM);
          [x1,x2] = ndgrid(1:DIM(1),1:DIM(2));
          for x3 = 1:DIM(3)
              M2  = inv(M\xM.VM(i).mat);
              y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*x3+M2(1,4));
              y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*x3+M2(2,4));
              y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*x3+M2(3,4));
              v(:,:,x3) = spm_bsplins(C, y1,y2,y3, [0 0 0 0 0 0]') > 0;
          end
          mask = mask & v;
          clear C v x1 x2 x3 M2 y1 y2 y3
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

    %-Cycle over bunches blocks within planes to avoid memory problems
    %==========================================================================
    swe_progress_bar('Init',nbchunks,'Parameter estimation','Chunks');

    for iChunk=1:nbchunks
      chunk = chunks(iChunk):chunks(iChunk+1)-1;
      
      %-Report progress
      %======================================================================
      if iChunk > 1, fprintf(repmat(sprintf('\b'),1,72)); end                  %-# 
      fprintf('%-40s: %30s', sprintf('Chunk %3d/%-3d',iChunk,nbchunks),...
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
      cmask(cmask) = any(diff(Y(:,cmask),1));        %-Mask constant data

      %-Mask out voxels where data is constant in at least one separable
      % matrix design either in a visit category or within-subject (BG - 27/05/2016)
      %------------------------------------------------------------------
      for g = 1:nGr_dof % first look data for each separable matrix design
        if sum(iGr_dof'==g) > 1 % do not look for cases where the separable matrix design is only one row (BG - 05/08/2016)     
          cmask(cmask)  = any(abs(diff(Y(iGr_dof'==g, cmask),1)) > eps, 1); % mask constant data within separable matrix design g (added by BG on 29/08/16)
          if isfield(SwE.type,'modified') % added by BG on 29/08/16
            for g2 = 1:nGr % then look data for each "homogeneous" group
              % check if the data is contant over subject for each visit category
              for k = 1:nVis_g(g2) 
                if sum(iGr_dof'==g & iGr == uGr(g2) & iVis == uVis_g{g2}(k)) > 1 % do not look for cases when the data is only one row (BG - 05/08/2016)
                  cmask(cmask) = any(abs(diff(Y(iGr_dof'==g & iGr == uGr(g2) & iVis == uVis_g{g2}(k), cmask),1)) > eps, 1);
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
      clear diffVis subjList

      Y = Y(:, cmask); %-Data within mask
      CrS = sum(cmask);

      %==================================================================
      %-Proceed with General Linear Model (if there are voxels)
      %==================================================================
      if CrS
        beta  = pX*Y;                     %-Parameter estimates
        if SwE.SS >= 4  % Cluster-wise adjustments
            res = zeros(size(Y));
            for i = 1:nSubj
                res(iSubj==uSubj(i),:) = corr{i} *...
                    (Y(iSubj==uSubj(i),:)-xX.X(iSubj==uSubj(i),:)*beta);
            end
        else
            res = diag(corr)*(Y-xX.X*beta); %-Corrected residuals
        end
        clear Y                           %-Clear to save memory

        %-Estimation of the data variance-covariance components (modified SwE) 
        %-SwE estimation (classic version)
        %--------------------------------------------------------------
        c = zeros(numel(chunk),1);

        if isfield(SwE.type,'modified')
            Cov_vis=zeros(nCov_vis,CrS);
            for i = Ind_Cov_vis_diag
                Cov_vis(i,:) = mean(res(Flagk(i,:),:).^2, 1);
            end
            % Check if some voxels have variance < eps and mask them 
            tmp = ~any(Cov_vis(Ind_Cov_vis_diag,:) < eps); % modified by BG on 29/08/16
            if any(~tmp)
                beta    = beta(:,tmp);
                res     = res(:,tmp);
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
                
                % compute the beta covariance matrice(s)
                switch dof_type 
                  case 1
                    Cov_beta = zeros(nCov_beta, CrS);
                    it = 0;
                    for g = 1:nGr
                        Cov_beta_g = weight(:,iGr_Cov_vis_g==g) * Cov_vis(iGr_Cov_vis_g==g,:);
                        for i=1:nCov_beta
                            it = it + 1;
                            c(cmask) = Cov_beta_g(i,:);
                            Vcov_beta_g(it) = swe_data_write(Vcov_beta_g(it), c, chunk);
                        end
                        Cov_beta = Cov_beta + Cov_beta_g;
                    end
                  case {0 2 3}
                    Cov_beta = weight * Cov_vis;
                end
            end
        else % else for "if isfield(SwE.type,'modified')"
            Cov_beta = 0;
            it = 0;
            for i = 1:nSubj
                Cov_beta_i_tmp = weight(:,Ind_Cov_vis_classic==i) *...
                    (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
                Cov_beta = Cov_beta + Cov_beta_i_tmp;
                if dof_type == 1 %need to save all subject contributions...
                  for ii=1:nCov_beta
                    it = it + 1;
                    c(cmask) = Cov_beta_i_tmp(ii,:);
                    Vcov_beta_g(it) = swe_data_write(Vcov_beta_g(it), c, chunk);
                  end
                end
            end
        end
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
      for iBeta=1:nBeta
        if CrS
          c(cmask) = beta(iBeta,:);
        end
        Vbeta(iBeta) = swe_data_write(Vbeta(iBeta), c, chunk); 
      end

      %-Write CovVis files if needed
      %----------------------------------------------------------------------
      if isfield(SwE.type,'modified') 
        for iCov_vis=1:nCov_vis
          if CrS
            c(cmask) = Cov_vis(iCov_vis,:);
          end
          Vcov_vis(iCov_vis) = swe_data_write(Vcov_vis(iCov_vis), c, chunk);
        end
      end

      %-Write CovBeta files
      %----------------------------------------------------------------------
      for iCov_beta=1:nCov_beta
        if CrS
          c(cmask) = Cov_beta(iCov_beta,:);
        end
        Vcov_beta(iCov_beta) = swe_data_write(Vcov_beta(iCov_beta), c, chunk);
      end
      
      %-Report progress
      %======================================================================
      fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');             %-#
      swe_progress_bar('Set',i);
    end % iChunk=1:nbchunks

    swe_progress_bar('Clear');        
    
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
    
else % matrix input
    % check how the data image treat 0 (as NaN or not)
    YNaNrep = 0;
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
    for g = 1:nGr_dof % first look data for each separable matrix design
      if sum(iGr_dof'==g) > 1 % do not look for cases where the separable matrix design is only one row (BG - 05/08/2016)     
        cmask(cmask) = any(abs(diff(Y(iGr_dof'==g,cmask),1)) > eps, 1); % mask constant data within separable matrix design g (added by BG on 29/08/16)
        if isfield(SwE.type,'modified') % added by BG on 29/08/16
          for g2 = 1:nGr % then look data for each "homogeneous" group
            % check if the data is contant over subject for each visit category
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
    clear diffVis

    Y      = Y(:,cmask);                          %-Data within mask
    CrS    = sum(cmask);                          %-# current voxels

    %==================================================================
    %-Proceed with General Linear Model (if there are voxels)
    %==================================================================
    if CrS

        %-General linear model: Ordinary least squares estimation
        %--------------------------------------------------------------
        fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...estimation');%-#

        crBeta  = pX*Y;                     %-Parameter estimates

        if SwE.SS >= 4  % Cluster-wise adjustments
            res = zeros(size(Y));
            for i = 1:nSubj
                res(iSubj==uSubj(i),:) = corr{i} *...
                    (Y(iSubj==uSubj(i),:)-xX.X(iSubj==uSubj(i),:)*crBeta);
            end
        else
            res   = diag(corr)*(Y-xX.X*crBeta); %-Corrected residuals
        end
        clear Y                        %-Clear to save memory

        %-Estimation of the data variance-covariance components (modified SwE)
        %-SwE estimation (classic version)
        %--------------------------------------------------------------

        if isfield(SwE.type,'modified')
            crCov_beta = 0;
            crCov_vis=zeros(nCov_vis,CrS);
            for i = Ind_Cov_vis_diag
                crCov_vis(i,:) = mean(res(Flagk(i,:),:).^2, 1);
            end
            % Check if some voxels have variance < eps and mask them
            tmp = ~any(crCov_vis(Ind_Cov_vis_diag,:) < eps); % modified by BG on 29/08/16
            if any(~tmp)
                crBeta    = crBeta(:,tmp);
                res     = res(:,tmp);
                cmask(cmask)  = tmp;
                CrS     = sum(cmask);
                crCov_vis = crCov_vis(:,tmp);
            end
            if CrS % Check if there is at least one voxel left
                for i = Ind_Cov_vis_off_diag
                    if any(Flagk(i,:))
                        crCov_vis(i,:)= sum(res(Flagk(i,:),:).*res(Flagkk(i,:),:), 1).*...
                            sqrt(crCov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,1)),:).*...
                            crCov_vis(Ind_Cov_vis_diag(Ind_corr_diag(i,2)),:)./...
                            sum(res(Flagk(i,:),:).^2, 1)./...
                            sum(res(Flagkk(i,:),:).^2, 1));
                    end
                end
                %NaN may be produced in crCov. estimation when one correspondant
                %variance are = 0, so set them to 0
                crCov_vis(isnan(crCov_vis))=0;
                %need to check if the eigenvalues of crCov_vis matrices are >=0
                for g = 1:nGr
                    for iVox = 1:CrS
                        tmp = zeros(nVis_g(g));
                        tmp(tril(ones(nVis_g(g)))==1) = crCov_vis(iGr_Cov_vis_g==g,iVox);
                        tmp = tmp + tmp' - diag(diag(tmp));
                        [V, D] = eig(tmp);
                        if any (diag(D)<0) %Bug corrected (BG - 19/09/13)
                            D(D<0) = 0;
                            tmp = V * D * V';
                            crCov_vis(iGr_Cov_vis_g==g,iVox) = tmp(tril(ones(nVis_g(g)))==1); %Bug corrected (BG - 19/09/13)
                        end
                    end
                end
            end
        else % else for "if isfield(SwE.type,'modified')"
            if dof_type == 1 %need to save all subject contributions...
                crCov_beta_i =  zeros(nSubj,nCov_beta,CrS);
            end
            crCov_beta = 0;
            for i = 1:nSubj
                crCov_beta_i_tmp = weight(:,Ind_Cov_vis_classic==i) *...
                    (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
                crCov_beta = crCov_beta + crCov_beta_i_tmp;
                if dof_type == 1 %need to save all subject contributions...
                    crCov_beta_i(i,:,:) = crCov_beta_i_tmp;
                end
            end
        end
        clear res
        if isfield(SwE.type,'modified')
            fprintf('\n');                                                    %-#
            disp('Working on the SwE computation...');
            switch dof_type 
                case 1
                    crCov_beta = zeros(nCov_beta,CrS); % initialize SwE for the plane
                    crCov_beta_i =  zeros(nGr,nCov_beta,CrS);
                    for g = 1:nGr
                        crCov_beta_i(g,:,:) = weight(:,iGr_Cov_vis_g==g) * crCov_vis(iGr_Cov_vis_g==g,:);                
                        Cov_beta = Cov_beta + crCov_beta_i(g,:,:);
                    end
                case {0 2 3}
                    crCov_beta = weight * crCov_vis;
            end
        end
        fprintf('\n');                                                    %-#
    end
        
    %-Save betas etc. 
    %----------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...saving results'); %-#

    mask = cmask;       
    save(sprintf('swe_%s_mask%s',file_data_type,file_ext), 'mask');
    clear mask

    beta = zeros(nBeta, nVox);
    if CrS
      beta(:,cmask) = crBeta;
    end
    save(sprintf('swe_%s_beta_b%s',file_data_type,file_ext), 'beta');
    clear beta crBeta

    if isfield(SwE.type,'modified')
        cov_vis = zeros(nCov_vis, nVox);
        if CrS
          cov_vis(:,cmask) = crCov_vis;
        end
        save(sprintf('swe_%s_cov_vv%s',file_data_type,file_ext), 'cov_vis');
        clear cov_vis crCov_vis
    end

    cov_beta = zeros(nCov_beta, nVox);
    cov_beta(:,cmask) = crCov_beta;
    save(sprintf('swe_%s_cov%s',file_data_type,file_ext), 'cov_beta');
    clear cov_beta crCov_beta
    if dof_type == 1
        nGr = nSubj;
        cov_beta_g = zeros(nGr, nCov_beta, nVox);
        cov_beta_g(:,:,cmask) = crCov_beta_i;
        save(sprintf('swe_%s_cov_g_bb%s',file_data_type,file_ext), 'cov_beta_g');
        clear cov_beta_g crCov_beta_i
    end
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done');   
    
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...some clean up'); %-#
  
    XYZ         = [];
    M           = [];
    DIM         = [];
    S           = CrS;
    VM          = sprintf('swe_%s_mask%s', file_data_type, file_ext);
    Vbeta       = sprintf('swe_%s_beta_b%s', file_data_type, file_ext);
    Vcov_beta   = sprintf('swe_%s_cov%s', file_data_type, file_ext);
    Vcov_vis    = sprintf('swe_%s_cov_vv%s', file_data_type, file_ext);
    Vcov_beta_g = sprintf('swe_%s_cov_g_bb%s', file_data_type, file_ext);
end

%-place fields in SwE
%--------------------------------------------------------------------------
SwE.xVol.XYZ   = XYZ;               %-InMask XYZ coords (voxels)
SwE.xVol.M     = M;                 %-voxels -> mm
SwE.xVol.iM    = inv(M);            %-mm -> voxels
SwE.xVol.DIM   = DIM';               %-image dimensions
% SwE.xVol.FWHM  = FWHM;              %-Smoothness data
% SwE.xVol.R     = R;                 %-Resel counts
SwE.xVol.S     = S;                 %-Volume (voxels)
% SwE.xVol.VRpv  = VRpv;              %-Filehandle - Resels per voxel
SwE.xVol.units = {'mm' 'mm' 'mm'};

SwE.Vbeta      = Vbeta;             %-Filehandle - Beta
SwE.Vcov_beta  = Vcov_beta;         %-Filehandle - Beta covariance
if isfield(SwE.type,'modified')
    SwE.Vcov_vis   = Vcov_vis;      %-Filehandle - Visit covariance    
end
if dof_type == 1
    SwE.Vcov_beta_g  = Vcov_beta_g;     %-Filehandle - Beta covariance contributions
end
% if ~all(vFWHM==0)
%     SwE.Vsmcov_vis = Vsmcov_vis;    %-Filehandle - Visit covariance
% end
SwE.VM         = VM;                %-Filehandle - Mask

SwE.xX         = xX;                %-design structure
SwE.xM         = xM;                %-mask structure

SwE.xCon       = struct([]);        %-contrast structure

SwE.swd        = pwd;

SwE.Subj.uSubj = uSubj;
SwE.Subj.nSubj = nSubj;

if isfield(SwE.type,'modified')
    
    if dof_type >1
        SwE.Vis.weight        = weight;
        SwE.Vis.iGr_Cov_vis_g = iGr_Cov_vis_g;
        SwE.Vis.Ind_corr_diag = Ind_corr_diag;
    end
    SwE.Vis.uVis_g = uVis_g;
    SwE.Vis.nVis_g = nVis_g;
    SwE.Vis.nCov_vis_g = nCov_vis_g;
    SwE.Vis.nCov_vis = nCov_vis;
    SwE.Vis.max_nVis_g = max_nVis_g;
    SwE.Vis.min_nVis_g = min_nVis_g;
    
    SwE.Gr.uGr       = uGr;
    SwE.Gr.nGr       = nGr;
    SwE.Gr.nSubj_g   = nSubj_g;
    SwE.Gr.uSubj_g   = uSubj_g;
else
    if dof_type == 1
        SwE.Gr.nGr   = nSubj;
    end
end

SwE.dof.uGr_dof   = uGr_dof; 
SwE.dof.nGr_dof   = nGr_dof;
SwE.dof.iGr_dof   = iGr_dof; 
SwE.dof.iBeta_dof = iBeta_dof;
SwE.dof.pB_dof    = pB_dof;
SwE.dof.nSubj_dof = nSubj_dof;
SwE.dof.edof_Subj = edof_Subj;
SwE.dof.dof_type  = dof_type;
if dof_type == 1 
    SwE.dof.edof_Gr = edof_Gr;
elseif dof_type == 0
    SwE.dof.dof_cov = dof_cov;
else
    SwE.dof.dofMat = dofMat; 
end

% save the version number of the toolbox
SwE.ver = swe('ver');

%-Save analysis parameters in SwE.mat file
%--------------------------------------------------------------------------
fprintf('%-40s: %30s','Saving SwE.mat','...writing');                   %-#
if isOctave
    save('SwE.mat','SwE');
elseif spm_matlab_version_chk('7') >=0
    save('SwE','SwE','-V6');
else
    save('SwE','SwE');
end
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#

%==========================================================================
%- E N D: Cleanup GUI
%==========================================================================
%spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
fprintf('...use the results section for assessment\n\n')  
