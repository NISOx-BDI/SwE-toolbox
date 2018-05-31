function swe_cp(SwE)
 
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
if exist(fullfile(SwE.swd,'swe_vox_mask.nii'),'file') == 2
 
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
 
files = {'^swe_vox_mask\..{3}$','^swe_vox_beta_bb\..{3}$','^swe_vox_cov\..{3}$',...
    '^swe_vox_cov_vv.{4}\..{3}$','^swe_vox_edf_c.{4}\..{3}$','^swe_vox_\w{1}stat_c.{4}\..{3}$',...
    '^swe_vox_\w{2}stat_c.{4}\..{3}$'};
 
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
        tmp(i,:) = any(xX.X(iGr_dof==i,:));
    end
    if nGr_dof==1 | all(sum(tmp,1)==1) %#ok<OR2>
        break % all is ok, just stop the while
    else
        ind1 = find(sum(tmp,1)>1,1); % detect the first column in common
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
    DIM      = VY(1).dim(1:3)';
    % VOX      = sqrt(diag(M(1:3, 1:3)'*M(1:3, 1:3)))'; (commented by BG on 08/11/2016)
    xdim     = DIM(1); ydim = DIM(2); zdim = DIM(3);
    %vFWHM    = SwE.vFWHM; to be added later (for the variance smoothing)

    % check how the data image treat 0 (as NaN or not)
    YNaNrep = VY(1).dt(2);

    %-Maximum number of residual images for smoothness estimation
    %--------------------------------------------------------------------------
    % MAXRES   = Inf; (commented by BG on 08/11/2016)
    % nSres    = nScan; (commented by BG on 08/11/2016)

    fprintf('%-40s: %30s','Output images','...initialising');           %-#

    %-Initialise new mask name: current mask & conditions on voxels
    %----------------------------------------------------------------------
    disp(file_ext)
    VM    = swe_create_vol(sprintf('swe_vox_mask%s', file_ext), DIM, M,...
                           'swe_cp:resultant analysis mask', isMeshData);

    %-Initialise beta image files
    %----------------------------------------------------------------------

    for i = 1:nBeta
        Vbeta(i) = swe_create_vol(sprintf('swe_vox_beta_b%02d%s',i,file_ext),...
                                  DIM, M,...
                                  sprintf('swe_cp:beta (%02d) - %s',i,xX.name{i}),...
                                  isMeshData);
    end

    %-Initialise Cov_beta image files
    %----------------------------------------------------------------------

    it=0;
    for i=1:nBeta
        for ii=i:nBeta
            it=it+1;
            Vcov_beta(it) = swe_create_vol(sprintf('swe_vox_cov_b%02d_b%02d%s',i,ii,file_ext),...
                                           DIM, M, sprintf('cov_beta_%02d_%02d hats - %s/%s',...
                                                i,ii,xX.name{i},xX.name{ii}),...
                                           isMeshData);
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
                    Vcov_beta_g(it) = swe_create_vol([sprintf('swe_vox_cov_g%02d_b%02d_b%02d',g,ii,iii) file_ext],...
                        DIM, M, sprintf('cov_beta_g_%02d_%02d_%02d hats - group %s - %s/%s',...
                            g,ii,iii,num2str(uGr(g)),xX.name{ii},xX.name{iii}), isMeshData);
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
                    Vcov_vis(it) = swe_create_vol([sprintf('swe_vox_cov_g%02d_v%02d_v%02d',g,ii,iii) file_ext],...
                                                  DIM, M, sprintf('cov_vis_%02d_%02d_%02d hats - group %s - visits %s/%s',...
                                                       g,ii,iii,num2str(uGr(g)),num2str(uVis_g{g}(ii)),num2str(uVis_g{g}(iii))),...
                                                       isMeshData);
                end
            end
        end
    end
    %-Initialise standardised residual images
    %----------------------------------------------------------------------
    % for i = 1:nSres
    %     VResI(i) = swe_create_vol(sprintf('swe_vox_resid_y%02d.img', i),...
    %                               DIM, M, sprintf('spm_spm:ResI (%02d)', i),...
    %                               isMeshData);
    % end
    % fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');    %-# 
    %  
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
    % use chunks instead of slices
    % chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
    % nbchunks  = ceil(prod(DIM) / chunksize);
    % chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

    %spm_progress_bar('Init',nbchunks,'Parameter estimation','Chunks');

    %-Initialise variables used in the loop
    %==========================================================================
    [xords, yords] = ndgrid(1:xdim, 1:ydim);
    xords = xords(:)'; yords = yords(:)';           % plane X,Y coordinates
    S     = 0;                                      % Volume (voxels)
    % i_res = round(linspace(1,nScan,nSres))';        % Indices for residual (commented by BG on 08/11/2016)

    %-Initialise XYZ matrix of in-mask voxel co-ordinates (real space)
    %--------------------------------------------------------------------------
    XYZ   = zeros(3,xdim*ydim*zdim);

    %-Cycle over bunches blocks within planes to avoid memory problems
    %==========================================================================
    str   = 'parameter estimation';
    spm_progress_bar('Init',100,str,'');

    for z = 1:nbz:zdim                       %-loop over planes (2D or 3D data)

        % current plane-specific parameters
        %----------------------------------------------------------------------
        CrPl         = z:min(z+nbz-1,zdim);       %-plane list
        zords        = CrPl(:)*ones(1,xdim*ydim); %-plane Z coordinates
        CrBl         = [];                        %-parameter estimates
    %     CrResI       = [];                        %-residuals (commented by BG on 08/11/2016)
        Q            = [];                        %-in mask indices for this plane
        if isfield(SwE.type,'modified')
            CrCov_vis    = []; 
        else
            CrCov_beta   = [];
            if dof_type ==1
                CrCov_beta_i = [];
            end
        end

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
                if SwE.SS >= 4  % Cluster-wise adjustments
                    res = zeros(size(Y));
                    for i = 1:nSubj
                        res(iSubj==uSubj(i),:) = corr{i} *...
                            (Y(iSubj==uSubj(i),:)-xX.X(iSubj==uSubj(i),:)*beta);
                    end
                else
                    res   = diag(corr)*(Y-xX.X*beta); %-Corrected residuals
                end
                clear Y                           %-Clear to save memory

                %-Estimation of the data variance-covariance components (modified SwE) 
                %-SwE estimation (classic version)
                %--------------------------------------------------------------

                if isfield(SwE.type,'modified')
                    Cov_beta = 0;
                    Cov_vis=zeros(nCov_vis,CrS);
                    for i = Ind_Cov_vis_diag
                        Cov_vis(i,:) = mean(res(Flagk(i,:),:).^2, 1);
                    end
                    % Check if some voxels have variance < eps and mask them 
                    tmp = ~any(Cov_vis(Ind_Cov_vis_diag,:) < eps); % modified by BG on 29/08/16
                    if any(~tmp)
                        beta    = beta(:,tmp);
                        res     = res(:,tmp);
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
                else % else for "if isfield(SwE.type,'modified')"
                    if dof_type == 1 %need to save all subject contributions...
                        Cov_beta_i =  NaN(nSubj,nCov_beta,CrS);
                    end
                    Cov_beta = 0;
                    for i = 1:nSubj
                        Cov_beta_i_tmp = weight(:,Ind_Cov_vis_classic==i) *...
                            (res(Indexk(Ind_Cov_vis_classic==i),:) .* res(Indexkk(Ind_Cov_vis_classic==i),:));
                        Cov_beta = Cov_beta + Cov_beta_i_tmp;
                        if dof_type == 1 %need to save all subject contributions...
                            Cov_beta_i(i,:,:) = Cov_beta_i_tmp;
                        end
                    end
                end

                %-Save betas etc. for current plane as we go along
                %----------------------------------------------------------
                CrBl          = [CrBl,    beta]; %#ok<AGROW>
    %             CrResI        = [CrResI,  res(i_res,:)]; %#ok<AGROW> (commented by BG on 08/11/2016)
                if isfield(SwE.type,'modified') 
                    CrCov_vis     = [CrCov_vis,  Cov_vis]; %#ok<AGROW>
                else
                    CrCov_beta     = [CrCov_beta, Cov_beta]; %#ok<AGROW>
                    if dof_type == 1
                        CrCov_beta_i     = cat(3, CrCov_beta_i, Cov_beta_i);
                    end
                end
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

        %-Write beta images
        %------------------------------------------------------------------
        for i = 1:nBeta
            if ~isempty(Q), jj(Q) = CrBl(i,:); end
            Vbeta(i) = spm_write_plane(Vbeta(i), jj, CrPl);
        end

        %-Write visit covariance images
        %------------------------------------------------------------------
        if isfield(SwE.type,'modified')
            for i=1:nCov_vis
                if ~isempty(Q), jj(Q) = CrCov_vis(i,:); end
                Vcov_vis(i) = spm_write_plane(Vcov_vis(i), jj, CrPl);
            end
        end

        %-Write SwE images and contributions if needed
        %------------------------------------------------------------------
        if isfield(SwE.type,'classic')       
            for i=1:nCov_beta
                if ~isempty(Q), jj(Q) = CrCov_beta(i,:); end
                Vcov_beta(i) = spm_write_plane(Vcov_beta(i), jj, CrPl);
            end
            if dof_type == 1
                it = 0;
                for i=1:nSubj
                    for ii=1:nCov_beta
                        it = it + 1;
                        if ~isempty(Q), jj(Q) = CrCov_beta_i(i,ii,:); end
                        Vcov_beta_g(it) = spm_write_plane(Vcov_beta_g(it), jj, CrPl);
                    end
                end
            end
        end

        %-Write standardised residual images
        %------------------------------------------------------------------
    %     for i = 1:nSres
    %         if ~isempty(Q), jj(Q) = CrResI(i,:)./...
    %                 sqrt(CrCov_vis(Flagk(:,i) & Flagkk(:,i),:)); 
    %         end 
    %         VResI(i) = spm_write_plane(VResI(i), jj, CrPl);
    %     end

        %-Report progress
        %----------------------------------------------------------------------
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done');   
        spm_progress_bar('Set',100*(bch + nbch*(z - 1))/(nbch*zdim));

    end % (for z = 1:zdim)
    fprintf('\n');                                                          %-#
    spm_progress_bar('Clear')
    clear beta res Cov_vis CrBl CrResI CrCov_vis jj%-Clear to save memory
    if isfield(SwE.type,'modified')
        clear Cov_vis CrCov_vis
    else
        clear  Cov_beta CrCov_beta        
        if dof_type == 1
            clear Cov_beta_i CrCov_beta_i
        end
    end
    XYZ   = XYZ(:,1:S); % remove all the data not used 

    %-SwE computation (for modified version, done later in case of a spatial regul.)
    %==========================================================================
    if isfield(SwE.type,'modified')

        %-Loading the visit covariance for the whole brain
        %----------------------------------------------------------------------
        fprintf('Loading the visit covariance for the SwE computation...'); %-#

        Cov_vis = spm_get_data(Vcov_vis,XYZ);

        %- Spatial regularization of the visit covariance if required
        %----------------------------------------------------------------------
        % Blurred mask is used to truncate kernel to brain; if not
        % used variance at edges would be underestimated due to
        % convolution with zero activity out side the brain.
        %-----------------------------------------------------------------
    %     Q           = cumprod([1,DIM(1:2)'])*XYZ - ...
    %         sum(cumprod(DIM(1:2)'));
    %     if ~all(vFWHM==0)
    %         fprintf('Working on the SwE spatial regularization...'); %-#
    %         SmCov_vis = zeros(xdim, ydim, zdim);
    %         SmMask    = zeros(xdim, ydim, zdim);
    %         TmpVol    = zeros(xdim, ydim, zdim);
    %         TmpVol(Q) = ones(size(Q));
    %         spm_smooth(TmpVol,SmMask,vFWHM./VOX);
    %         jj = NaN(xdim,ydim,zdim);
    %         for i = 1:nCov_vis
    %             TmpVol(Q) = Cov_vis(i,:);
    %             spm_smooth(TmpVol,SmCov_vis,vFWHM./VOX);
    %             Cov_vis (i,:) = SmCov_vis(Q)./SmMask(Q);
    %             jj(Q) = Cov_vis (i,:);
    %             spm_write_vol(Vsmcov_vis(i),jj);
    %         end
    %     end
        fprintf('\n');                                                    %-#
        disp('Working on the SwE computation...');
        %Computation of the SwE
        str   = 'SwE computation';
        spm_progress_bar('Init',100,str,'');

        S_z = 0;
        for z = 1:zdim                       %-loop over planes (2D or 3D data)       
            XY_z = XYZ(1:2,XYZ(3,:)==z); % extract coord in plane z        
            Q_z = cumprod([1,DIM(1)'])*XY_z - ...
                sum(cumprod(DIM(1)'));
            s_z = length(Q_z); % number of active voxels in plane z
            jj = NaN(xdim,ydim);
            switch dof_type 
                case 1
                    Cov_beta = zeros(nCov_beta,s_z); % initialize SwE for the plane
                    it = 0;
                    for g = 1:nGr
                        Cov_beta_g = weight(:,iGr_Cov_vis_g==g) * Cov_vis(iGr_Cov_vis_g==g,(1+S_z):(S_z+s_z));
                        for i=1:nCov_beta
                            if ~isempty(Q_z), jj(Q_z)=Cov_beta_g(i,:); end
                            it = it + 1;
                            Vcov_beta_g(it)=spm_write_plane(Vcov_beta_g(it),jj, z);
                        end
                        Cov_beta = Cov_beta + Cov_beta_g;
                        spm_progress_bar('Set',100*((z-1)/zdim + g/nGr/zdim));
                    end
                case {0 2 3}
                    Cov_beta = weight * Cov_vis(:,(1+S_z):(S_z+s_z));
                    spm_progress_bar('Set',100*(z/zdim));
            end
            for i=1:nCov_beta
                if ~isempty(Q_z), jj(Q_z)=Cov_beta(i,:); end
                Vcov_beta(i)=spm_write_plane(Vcov_beta(i),jj, z);
            end
            S_z = S_z + s_z;
        end% (for z = 1:zdim)
        fprintf('\n');                                                    %-#
        spm_progress_bar('Clear')

    end


    %==========================================================================
    % - P O S T   E S T I M A T I O N   C L E A N U P
    %==========================================================================
    if S == 0, spm('alert!','No inmask voxels - empty analysis!'); return; end

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

    Y      = Y(:,Cm);                          %-Data within mask
    CrS    = sum(Cm);                          %-# current voxels

    %==================================================================
    %-Proceed with General Linear Model (if there are voxels)
    %==================================================================
    if CrS

        %-General linear model: Ordinary least squares estimation
        %--------------------------------------------------------------
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...estimation');%-#

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
                Cm(Cm)  = tmp;
                CrS     = sum(Cm);
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
                crCov_beta_i =  NaN(nSubj,nCov_beta,CrS);
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
                    crCov_beta_i =  NaN(nGr,nCov_beta,CrS);
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

    mask = Cm;       
    save(sprintf('swe_vox_mask%s',file_ext), 'mask');
    clear mask

    beta = NaN(nBeta, nVox);
    beta(:,Cm) = crBeta;
    save(sprintf('swe_vox_beta_bb%s',file_ext), 'beta');
    clear beta crBeta

    if isfield(SwE.type,'modified')
        cov_vis = NaN(nCov_vis, nVox);
        cov_vis(:,Cm) = crCov_vis;
        save(sprintf('swe_vox_cov_vv%s',file_ext), 'cov_vis');
        clear cov_vis crCov_vis
    end

    cov_beta = NaN(nCov_beta, nVox);
    cov_beta(:,Cm) = crCov_beta;
    save(sprintf('swe_vox_cov%s',file_ext), 'cov_beta');
    clear cov_beta crCov_beta
    if dof_type == 1
        nGr = nSubj;
        cov_beta_g = NaN(nGr, nCov_beta, nVox);
        cov_beta_g(:,:,Cm) = crCov_beta_i;
        save(sprintf('swe_vox_cov_g_bb%s',file_ext), 'cov_beta_g');
        clear cov_beta_g crCov_beta_i
    end
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done');   
    
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...some clean up'); %-#
  
    XYZ         = [];
    M           = [];
    DIM         = [];
    S           = CrS;
    VM          = sprintf('swe_vox_mask%s', file_ext);
    Vbeta       = sprintf('swe_vox_beta_bb%s', file_ext);
    Vcov_beta   = sprintf('swe_vox_cov%s', file_ext);
    Vcov_vis    = sprintf('swe_vox_cov_vv%s', file_ext);
    Vcov_beta_g = sprintf('swe_vox_cov_g_bb%s', file_ext);
end

%-place fields in SwE
%--------------------------------------------------------------------------
SwE.xVol.XYZ   = XYZ;               %-InMask XYZ coords (voxels)
SwE.xVol.M     = M;                 %-voxels -> mm
SwE.xVol.iM    = inv(M);            %-mm -> voxels
SwE.xVol.DIM   = DIM;               %-image dimensions
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
%- E N D: Cleanup GUI
%==========================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#
%spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
fprintf('...use the results section for assessment\n\n')  
