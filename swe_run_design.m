function out = swe_run_design(varargin)
%
% Sandwich Estimator Method for Neuroimaging Longitudinal Data Analysis, SwE.
% takes a harvested job data structure and rearranges data into SwE
% data structure, then saves SwE.mat file.
%
% INPUT
%   job    - harvested job data structure (see matlabbatch help)
%
% OUTPUT
%   out    - filename of saved data structure.

% Written by Bryan Guillaume

% Job variable
% -------------------------------------------------------------------------
job   = varargin{1};

%-Change directory
%--------------------------------------------------------------------------
original_dir = pwd;
cd(job.dir{1});

%-Ask about overwriting files from previous analyses...
%--------------------------------------------------------------------------
if exist(fullfile(job.dir{1},'SwE.mat'),'file')
    str = { 'Current directory contains existing SwE file:',...
        'Continuing will overwrite existing file!'};
    if spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
        fprintf('%-40s: %30s\n\n',...
            'Abort...   (existing SwE file)',spm('time'));
        return
    end
end

% If we've gotten to this point we're committed to overwriting files.
% Delete them so we don't get stuck 
%--------------------------------------------------------------------------
files = {'^mask\..{3}$','^ResMS\..{3}$','^RPV\..{3}$',...
    '^beta_.{4}\..{3}$','^con_.{4}\..{3}$','^ResI_.{4}\..{3}$',...
    '^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$',...
    '^cov_beta_.{4}_.{4}\..{3}$', '^cov_vis_.{4}_.{4}_.{4}\..{3}$',...
    '^edf_.{4}\..{3}$'};

for i=1:length(files)
    j = spm_select('List',pwd,files{i});
    for k=1:size(j,1)
        spm_unlink(deblank(j(k,:)));
    end
end

%-Option definitions
%==========================================================================

% %-Generic factor names
% %--------------------------------------------------------------------------
sF = {'sF1','sF2','sF3','sF4'};

%-Covariate by factor interaction options
sCFI = {'<none>';...                                        %-1
    'with sF1';'with sF2';'with sF3';'with sF4';...         %-2:5
    'with sF2 (within sF4)';'with sF3 (within sF4)'};       %-6,7

%-DesMtx argument components for covariate by factor interaction options
% (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
%--------------------------------------------------------------------------
CFIforms = {'[]',   'C',    '{}';...                        %-1
    'I(:,1)',       'FxC',  '{sF{1}}';...                   %-2
    'I(:,2)',       'FxC',  '{sF{2}}';...                   %-3
    'I(:,3)',       'FxC',  '{sF{3}}';...                   %-4
    'I(:,4)',       'FxC',  '{sF{4}}';...                   %-5
    'I(:,[4,2])',   'FxC',  '{sF{4},sF{2}}';...             %-6
    'I(:,[4,3])',   'FxC',  '{sF{4},sF{3}}' };              %-7

%-Centre (mean correction) options for covariates & globals            (CC)
% (options 9-12 are for centering of global when using AnCova GloNorm) (GC)
%--------------------------------------------------------------------------
sCC = {'around overall mean';...                            %-1
    'around sF1 means';...                                  %-2
    'around sF2 means';...                                  %-3
    'around sF3 means';...                                  %-4
    'around sF4 means';...                                  %-5
    'around sF2 (within sF4) means';...                     %-6
    'around sF3 (within sF4) means';...                     %-7
    '<no centering>';...                                    %-8
    'around user specified value';...                       %-9
    '(as implied by AnCova)';...                            %-10
    'GM';...                                                %-11
    '(redundant: not doing AnCova)'}';                      %-12

% %-DesMtx I forms for covariate centering options
% %--------------------------------------------------------------------------
% CCforms = {'ones(nScan,1)',CFIforms{2:end,1},''}';

%-Global calculation options                                       (GXcalc)
%--------------------------------------------------------------------------
sGXcalc  = {'omit';...                                      %-1
    'user specified';...                                    %-2
    'mean voxel value (within per image fullmean/8 mask)'}; %-3

%-Global normalization options                                    (GloNorm)
%--------------------------------------------------------------------------
sGloNorm = {'AnCova';...                                    %-1
    'AnCova by sF1';...                                     %-2
    'AnCova by sF2';...                                     %-3
    'AnCova by sF3';...                                     %-4
    'AnCova by sF4';...                                     %-5
    'AnCova by sF2 (within sF4)';...                        %-6
    'AnCova by sF3 (within sF4)';...                        %-7
    'proportional scaling';...                              %-8
    '<no global normalisation>'};                           %-9

%-Grand mean scaling options                                        (GMsca)
% (NB: Grand mean scaling by subject is redundent for proportional scaling)
%--------------------------------------------------------------------------
sGMsca = {'scaling of overall grand mean';...               %-1
    'scaling of sF1 grand means';...                        %-2
    'scaling of sF2 grand means';...                        %-3
    'scaling of sF3 grand means';...                        %-4
    'scaling of sF4 grand means';...                        %-5
    'scaling of sF2 (within sF4) grand means';...           %-6
    'scaling of sF3 (within sF4) grand means';...           %-7
    '(implicit in PropSca global normalisation)';...        %-8
    '<no grand Mean scaling>'   };                          %-9

%-Conditions of no interest defaults
%--------------------------------------------------------------------------
B      = [];
Bnames = {};
factor = [];
DesName = 'SwE';

P = job.scans;

n = length(job.subjects);
% check length of variables
if ~(length(P) == 1 || length(P) == n)
    error('The number of scans and the length of the subject indicator variable does not match.')
end
I = (1:n)';
I = [I,ones(n,3)];

factor(1).name     = '';
factor(1).levels   = 1;
factor(1).variance = 0;
factor(1).dept     = 0;
% Set up subjects information: Subj
%--------------------------------------------------------------------------

Subj.iSubj = job.subjects;              %subjects list
Subj.nSubj = length(unique(Subj.iSubj));  %number of subjects

% Set up visits & groups information: Vis & Gr
%--------------------------------------------------------------------------

switch char(fieldnames(job.type))

    case 'modified'
     
        Vis.iVis = job.type.modified.visits;
        Vis.nVis = length(unique(Vis.iVis));
        Gr.iGr   = job.type.modified.groups;
        Gr.nGr   = length(unique(Gr.iGr));
        SS       = job.type.modified.ss;
        if length(Vis.iVis) ~= n
            error('The lengths of the subject and visit indicator variables do not match.')
        end
        if length(Gr.iGr) ~= n
            error('The lengths of the subject and group indicator variables do not match.')
        end
    case 'classic'

        Vis =[];
        Gr =[];
        SS       = job.type.classic.ss;

end

nScan = size(I,1);

%-Covariate partition(s): interest (C) & nuisance (G) excluding global
%==========================================================================
dstr   = {'covariate','nuisance variable'};
C  = []; Cnames = [];  %-Covariate DesMtx partitions & names
G  = []; Gnames = [];
H  = []; Hnames = [];
B  = []; Bnames = [];
xC = [];                         %-Struct array to hold raw covariates

% Covariate options:
nc=length(job.cov); % number of covariates
for i=1:nc

    c      = job.cov(i).c;
    cname  = job.cov(i).cname;
    rc     = c;                         %-Save covariate value
    rcname = cname;                     %-Save covariate name
%     if job.cov(i).iCFI==1, % no interaction
%         iCFI=1;
%     else
%         % SPMs internal factor numbers are 1 higher than specified in user
%         % interface as, internally, the first factor is always `replication'
%         iCFI=job.cov(i).iCFI+1;
%     end
%     switch job.cov(i).iCC,
%         case 1
%             iCC=1;
%         case {2,3,4}
%             iCC=job.cov(i).iCC+1;
%         otherwise
%             iCC=job.cov(i).iCC+3;
%     end

    %-Centre within factor levels as appropriate
%     if any(iCC == [1:7]),
%         c = c - spm_meanby(c,eval(CCforms{iCC}));
%     end

    %-Do any interaction (only for single covariate vectors)
    %----------------------------------------------------------------------
%     if iCFI > 1             %-(NB:iCFI=1 if size(c,2)>1)
%         tI        = [eval(CFIforms{iCFI,1}),c];
%         tConst    = CFIforms{iCFI,2};
%         tFnames   = [eval(CFIforms{iCFI,3}),{cname}];
%         [c,cname] = spm_DesMtx(tI,tConst,tFnames);
%     elseif size(c,2)>1          %-Design matrix block
%         [null,cname] = spm_DesMtx(c,'X',cname);
%     else
         cname = {cname};
%     end

    %-Store raw covariate details in xC struct for reference
    %-Pack c into appropriate DesMtx partition
    %----------------------------------------------------------------------
    %-Construct description string for covariate
    str = {sprintf('%s',rcname)};
    if size(rc,2)>1, str = {sprintf('%s (block of %d covariates)',...
            str{:},size(rc,2))}; end
%     if iCC < 8, str=[str;{['used centered ',sCC{iCC}]}]; end
%     if iCFI> 1, str=[str;{['fitted as interaction ',sCFI{iCFI}]}]; end
    typ = 1;
    tmp = struct(...
        'rc',   rc,    'rcname', rcname,...
        'c',    c,     'cname',  {cname},...%         'iCC',  iCC,   'iCFI',   iCFI,...
        'type', typ,...
        'cols', [1:size(c,2)] + size([H,C],2) + size([B,G],2)*min(typ-1,1),...
        'descrip', {str});
    if isempty(xC), xC = tmp; else xC = [xC,tmp]; end

    C     = [C,c];

    Cnames = [Cnames; cname];

end
clear c tI tConst tFnames

%==========================================================================
% - C O N F I G U R E   D E S I G N -
%==========================================================================

%-Images & image info: Map Y image files and check consistency of
% dimensions and orientation / voxel size
%==========================================================================

fprintf('%-40s: ','Mapping files')    
P = job.scans;
[~,~,file_ext] = fileparts(P{1});
isMat = strcmpi(file_ext,'.mat');
if isMat
    VY = {};
else
    VY = spm_vol(char(P));
end
%-Check compatibility of images
%--------------------------------------------------------------------------
spm_check_orientations(VY);

fprintf('%30s\n','...done')  


%-Global values, scaling and global normalisation
%==========================================================================
%-Compute global values
%--------------------------------------------------------------------------

switch char(fieldnames(job.globalc))
    case 'g_omit',
        iGXcalc = 1;
    case 'g_user',
        iGXcalc = 2;
    case 'g_mean',
        iGXcalc = 3;
end

switch job.globalm.glonorm
    case 1,
        iGloNorm = 9;
    case 2,
        iGloNorm = 8;
    case 3,
        iGloNorm = 1;
end

factor(1).levels
if factor(1).levels > 1
    % Override if factor-specific ANCOVA has been specified
    for i=1:length(factor)
        if factor(i).ancova
            iGloNorm=i+2;
        end
    end
end

%-Analysis threshold mask
%--------------------------------------------------------------------------
%-Work out available options:
% -Inf=>None, real=>absolute, complex=>proportional, (i.e. times global)
M_T = -Inf;

switch char(fieldnames(job.masking.tm)),
    case 'tma',
        % Absolute
        M_T = job.masking.tm.tma.athresh;
    case 'tmr',
        % Relative
        M_T = job.masking.tm.tmr.rthresh*sqrt(-1);
        % Need to force calculation of globals
        if iGXcalc~=2, iGXcalc=3; end
    case 'tm_none'
        % None
        M_T = -Inf;
end

if iGXcalc==1 && (any(iGloNorm == [1:5 8]) || ...
        (factor(1).levels > 1 && any([factor.gmsca])))
    % Over-ride omission of global calculation if we need it
    disp(' ');
    disp('SPM needs estimates of global activity.');
    disp('But you have specified to omit this computation.');
    disp('SPM has overridden this omission and will automatically compute ');
    disp('globals as the mean value of within brain voxels.');
    disp(' ');
    iGXcalc = 3;
end
sGXcalc = sGXcalc{iGXcalc};

switch iGXcalc,
    case 1
        %-Don't compute => no GMsca (iGMsca==9) or GloNorm (iGloNorm==9)
        g = [];
    case 2
        %-User specified globals
        g = job.globalc.g_user.global_uval;
    case 3
        %-Compute as mean voxel value (within per image fullmean/8 mask)
        g = zeros(nScan,1);
        fprintf('%-40s: %30s','Calculating globals',' ')                %-#
        for i = 1:nScan
            str = sprintf('%3d/%-3d',i,nScan);
            fprintf('%s%30s',repmat(sprintf('\b'),1,30),str)            %-#
            g(i) = spm_global(VY(i)); % FIXME % for meshes
        end
        fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')        %-#
    otherwise
        error('illegal iGXcalc')
end
rg = g;

fprintf('%-40s: ','Design configuration')                               %-#

%-Grand mean scaling options                                        (GMsca)
%--------------------------------------------------------------------------
if iGloNorm==8
    iGMsca=8;   %-grand mean scaling implicit in PropSca GloNorm
else
    switch char(fieldnames(job.globalm.gmsca))
        case 'gmsca_yes',
            iGMsca=1;
        case 'gmsca_no',
            iGMsca=9;
    end
    if factor(1).levels > 1
        % Over-ride if factor-specific scaling has been specified
        for i=1:numel(factor)
            if factor(i).gmsca
                iGMsca=i+2;
            end
        end
    end
end

%-Value for PropSca / GMsca                                            (GM)
%--------------------------------------------------------------------------
switch iGMsca,
    case 9                                %-Not scaling (GMsca or PropSca)
        GM = 0;                           %-Set GM to zero when not scaling
    case 1                                %-Ask user value of GM
        GM = job.globalm.gmsca.gmsca_yes.gmscv;
    otherwise
        if iGloNorm==8
            switch char(fieldnames(job.globalm.gmsca))
                case 'gmsca_yes',
                    % Proportionally scale to this value
                    GM = job.globalm.gmsca.gmsca_yes.gmscv;
                case 'gmsca_no',
                    GM = 50;
            end
        else
            % Grand mean scaling by factor eg. scans are scaled so that the
            % mean global value over each level of the factor is set to GM
            GM=50;
        end
end

%-If GM is zero then don't GMsca! or PropSca GloNorm
if GM==0,
    iGMsca=9;
    if iGloNorm==8,
        iGloNorm=9;
    end
end

%-Sort out description strings for GloNorm and GMsca
%--------------------------------------------------------------------------
sGloNorm = sGloNorm{iGloNorm};
sGMsca   = sGMsca{iGMsca};
if iGloNorm==8
    sGloNorm = sprintf('%s to %-4g',sGloNorm,GM);
elseif iGMsca<8
    sGMsca   = sprintf('%s to %-4g',sGMsca,GM);
end

%-Scaling: compute global scaling factors gSF required to implement
% proportional scaling global normalisation (PropSca) or grand mean
% scaling (GMsca), as specified by iGMsca (& iGloNorm)
%--------------------------------------------------------------------------
switch iGMsca,
    case 8
        %-Proportional scaling global normalisation
        if iGloNorm~=8, error('iGloNorm-iGMsca(8) mismatch for PropSca'), end
        gSF    = GM./g;
        g      = GM*ones(nScan,1);
    case 1
        %-Grand mean scaling according to iGMsca
        gSF    = GM./spm_meanby(g,ones(nScan,1));
        g      = g.*gSF;
    case 9
        %-No grand mean scaling
        gSF    = ones(nScan,1);
    otherwise
        error('illegal iGMsca')
end

%-Apply gSF to memory-mapped scalefactors to implement scaling
%--------------------------------------------------------------------------
if ~strcmpi(file_ext,'.mat')
    for i = 1:nScan
        VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*gSF(i); % FIXME % for meshes
    end
end

%-Global centering (for AnCova GloNorm)                                (GC)
%-If not doing AnCova then GC is irrelevant
%--------------------------------------------------------------------------
if ~any(iGloNorm == [1:7])
    iGC = 12;
    gc  = [];
else
    iGC = 10;
    gc = 0;
end
%-AnCova: Construct global nuisance covariates partition (if AnCova)
%--------------------------------------------------------------------------
if any(iGloNorm == [1:7])

    %-Centre global covariate as requested
    %----------------------------------------------------------------------
    switch iGC, case {1,2,3,4,5,6,7}    %-Standard sCC options
        gc = spm_meanby(g,eval(CCforms{iGC}));
        case 8                  %-No centering
            gc = 0;
        case 9                  %-User specified centre
            %-gc set above
        case 10                 %-As implied by AnCova option
            gcB = spm_meanby(g,ones(nScan,1));
            gcW = spm_meanby(g,Subj.iSubj);
        case 11                 %-Around GM
            gc = GM;
        otherwise               %-unknown iGC
            error('unexpected iGC value')
    end

    %-AnCova - add scaled centred global to DesMtx `G' partition
    %----------------------------------------------------------------------
    fB = gcW - gcB;
    gnamesB = 'between-subject global';
    fW = g - gcW;
    gnamesW = 'within-subject global';

    %-Save GX info in xC struct for reference
    %----------------------------------------------------------------------
    str     = {sprintf('%s: %s',dstr{2},gnamesB)};
    if any(iGMsca==[1:7]), str=[str;{['(after ',sGMsca,')']}]; end
    if iGC ~= 8, str=[str;{['used centered ',sCC{iGC}]}]; end
    if iGloNorm > 1
        str=[str;{['fitted as interaction ',sCFI{iGloNorm}]}];
    end
    tmpB  = struct(  'rc',rg.*gSF,       'rcname',gnamesB,...
        'c',fB,          'cname' ,{gnamesB},...%'iCC',iGC,      'iCFI'  ,iGloNorm,...
        'type',         3,...
        'cols',1 + size([H C B G],2),...
        'descrip',      {str}       );
    str     = {sprintf('%s: %s',dstr{2},gnamesB)};
    if any(iGMsca==[1:7]), str=[str;{['(after ',sGMsca,')']}]; end
    if iGC ~= 8, str=[str;{['used centered ',sCC{iGC}]}]; end
    if iGloNorm > 1
        str=[str;{['fitted as interaction ',sCFI{iGloNorm}]}];
    end
    tmpW  = struct(  'rc',rg.*gSF,       'rcname',gnamesW,...
        'c',fW,          'cname' ,{gnamesW},...%'iCC',iGC,      'iCFI'  ,iGloNorm,...
        'type',         3,...
        'cols',2 + size([H C B G],2),...
        'descrip',      {str}       );
    
    G = [G,fB,fW]; Gnames = [Gnames; {gnamesB}; {gnamesW}];
    if isempty(xC), xC = [tmpB,tmpW]; else xC = [xC,tmpB,tmpW]; end

elseif iGloNorm==8 || iGXcalc>1

    %-Globals calculated, but not AnCova: Make a note of globals
    %----------------------------------------------------------------------
    if iGloNorm==8
        str = { 'global values: (used for proportional scaling)';...
            '("raw" unscaled globals shown)'};
    elseif isfinite(M_T) && ~isreal(M_T)
        str = { 'global values: (used to compute analysis threshold)'};
    else
        str = { 'global values: (computed but not used)'};
    end

    rcname ='global';
    tmp     = struct('rc',rg,    'rcname',rcname,...
        'c',{[]},   'cname' ,{{}},...%        'iCC',0,    'iCFI'  ,0,...
        'type',     3,...
        'cols',     {[]},...
        'descrip',  {str}           );

    if isempty(xC), xC = tmp; else xC = [xC,tmp]; end
end

%-Save info on global calculation in xGX structure
%--------------------------------------------------------------------------
xGX = struct(...
    'iGXcalc', iGXcalc,  'sGXcalc', sGXcalc,  'rg',rg,...
    'iGMsca',  iGMsca,   'sGMsca',  sGMsca,   'GM',GM,    'gSF',gSF,...
    'iGC',     iGC,      'sGC',     sCC{iGC}, 'gc',gc,...
    'iGloNorm',iGloNorm, 'sGloNorm',sGloNorm);

%-Make a description string
%--------------------------------------------------------------------------
if isinf(M_T)
    xsM.Analysis_threshold = 'None (-Inf)';
elseif isreal(M_T)
    xsM.Analysis_threshold = sprintf('images thresholded at %6g',M_T);
else
    xsM.Analysis_threshold = sprintf(['images thresholded at %6g ',...
        'times global'],imag(M_T));
end

%-Construct masking information structure and compute actual analysis
% threshold using scaled globals (rg.*gSF)
%--------------------------------------------------------------------------
if isreal(M_T),
    M_TH = M_T  * ones(nScan,1);    %-NB: -Inf is real
else
    M_TH = imag(M_T) * (rg.*gSF);
end

%-Implicit masking: Ignore zero voxels in low data-types?
%--------------------------------------------------------------------------
% (Implicit mask is NaN in higher data-types.)
if strcmpi(file_ext,'.mat')
    type = 16; % assume that there is a nan representation
else
    type = getfield(spm_vol(P{1,1}),'dt')*[1,0]';
end
if ~spm_type(type,'nanrep')
    M_I = job.masking.im;  % Implicit mask ?
    if M_I
        xsM.Implicit_masking = 'Yes: zero''s treated as missing';
    else
        xsM.Implicit_masking = 'No';
    end
else
    M_I = 1;
    xsM.Implicit_masking = 'Yes: NaN''s treated as missing';
end
%-Explicit masking
%--------------------------------------------------------------------------
if isempty(job.masking.em{:})
    VM = [];
    xsM.Explicit_masking = 'No';
else
    if isMat
        VM = job.masking.em;
        [~,~,file_ext_mask] = fileparts(VM{1});
        if ~strcmpi(file_ext_mask,'.mat')
            error('The explicit mask is not in ".mat" format as expected.')
        end
    else
        VM = spm_vol(char(job.masking.em));
    end
    xsM.Explicit_masking = 'Yes';
end

xM     = struct('T',M_T, 'TH',M_TH, 'I',M_I, 'VM',{VM}, 'xs',xsM);

%-Construct full design matrix (X), parameter names and structure (xX)
%==========================================================================
X      = [H C B G];
tmp    = cumsum([size(H,2), size(C,2), size(B,2), size(G,2)]);
xX     = struct(    'X',        X,...
    'iH',       [1:size(H,2)],...
    'iC',       [1:size(C,2)] + tmp(1),...
    'iB',       [1:size(B,2)] + tmp(2),...
    'iG',       [1:size(G,2)] + tmp(3),...
    'name',     {[Hnames; Cnames; Bnames; Gnames]},...
    'I',        I,...
    'sF',       {sF});


%-Design description (an nx2 cellstr) - for saving and display
%==========================================================================
tmp = {sprintf('%d condition, +%d covariate, +%d block, +%d nuisance',...
    size(H,2),size(C,2),size(B,2),size(G,2));...
    sprintf('%d total, having %d degrees of freedom',...
    size(X,2),rank(X));...
    sprintf('leaving %d degrees of freedom from %d images',...
    size(X,1)-rank(X),size(X,1))};
xsDes = struct('Design',    {DesName},...
    'Global_calculation',   {sGXcalc},...
    'Grand_mean_scaling',   {sGMsca},...
    'Global_normalisation', {sGloNorm},...
    'Parameters',           {tmp});


fprintf('%30s\n','...done')                                             %-#


%==========================================================================
% - WB configuration - Only if needed
%==========================================================================
if isfield(job.WB, 'WB_yes')
  
  WB.RWB            = job.WB.WB_yes.WB_type;
  WB.SS             = job.WB.WB_yes.WB_ss;
  WB.nB             = job.WB.WB_yes.WB_nB;
  WB.RSwE           = job.WB.WB_yes.WB_SwE;
  WB.voxelWise      = 1;
  WB.voxelWiseInfo = [];
  switch char(fieldnames(job.WB.WB_yes.WB_cluster))
    
    case 'WB_cluster_no'
      WB.clusterWise  = 0;
      
    case 'WB_cluster_yes'
      WB.clusterWise  = 1;
      WB.clusterInfo = [];
      WB.clusterInfo.primaryThreshold = job.WB.WB_yes.WB_cluster.WB_cluster_yes;
      if WB.clusterInfo.primaryThreshold > 1 || WB.clusterInfo.primaryThreshold < 0
        error('cluster-forming threshold should be between 0 an 1 (this is a probability)');
      end
      
    case 'WB_cluster_yes_mat'
      WB.clusterWise  = 1;
      WB.clusterInfo = [];
      WB.clusterInfo.primaryThreshold = job.WB.WB_yes.WB_cluster.WB_cluster_yes_mat.WB_cluster_yes_mat_clusP;
      if WB.clusterInfo.primaryThreshold > 1 || WB.clusterInfo.primaryThreshold < 0
        error('cluster-forming threshold should be between 0 an 1 (this is a probability)');
      end
      if job.WB.WB_yes.WB_cluster.WB_cluster_yes_mat.WB_cluster_yes_mat_type == 0
        WB.clusterInfo.Vxyz = job.WB.WB_yes.WB_cluster.WB_cluster_yes_mat.WB_cluster_yes_mat_loc;
      elseif job.WB.WB_yes.WB_cluster.WB_cluster_yes_mat.WB_cluster_yes_mat_type == 1
        WB.clusterInfo.Vfaces = job.WB.WB_yes.WB_cluster.WB_cluster_yes_mat.WB_cluster_yes_mat_loc;
      else
        error('wrong specification of the type of ".mat" files')
      end     
  end
  
  switch char(fieldnames(job.WB.WB_yes.WB_stat))
    
    case 'WB_T'
      WB.stat = 'T';
      WB.con = job.WB.WB_yes.WB_stat.WB_T.WB_T_con;
      if any(size(WB.con) ~= [1 size(X,2)])
        error('contrast not well specified');
      end
      
    case 'WB_F'
      WB.stat = 'F';
      WB.con = job.WB.WB_yes.WB_stat.WB_F.WB_F_con;
      if size(WB.con,2) ~= size(X,2)
        error('contrast not well specified');
      end
      
    otherwise
      error('unexpected statistic type');
  end
 
end


%-Assemble SwE structure like it is done in SPM structure
%==========================================================================
SwE.xY.P      = P;            % filenames
SwE.xY.VY     = VY;           % mapped data
SwE.xY.isMat  = isMat;
SwE.nscan     = size(xX.X,1); % scan number
SwE.xX        = xX;           % design structure
SwE.xC        = xC;           % covariate structure
SwE.xGX       = xGX;          % global structure
SwE.xM        = xM;           % mask structure
SwE.xsDes     = xsDes;        % description
SwE.type      = job.type;     % SwE type (modified or classic)
SwE.SS        = SS;           % SwE small samples adj. type
SwE.Subj      = Subj;         % subjects data
SwE.Vis       = Vis;          % visits data (empty if classic SwE)
SwE.Gr        = Gr;           % groups data (empty if classic SwE)
if isfield(job.WB, 'WB_yes')
  SwE.WB      = WB;           % WB structure
end

%-Save SwE.mat and set output argument
%--------------------------------------------------------------------------
fprintf('%-40s: ','Saving SwE configuration')                           %-#

if spm_check_version('matlab','7') >= 0
    save('SwE.mat', 'SwE', '-V6');
else
    save('SwE.mat', 'SwE');
end
fprintf('%30s\n','...SwE.mat saved')                                    %-#

out.swemat{1} = fullfile(pwd, 'SwE.mat');

%-Display Design report
%==========================================================================
if ~spm('CmdLine')
    fprintf('%-40s: ','Design reporting') 
    if strcmpi(file_ext,'.mat')
        fname = cellstr(repmat('  ', nScan, 1));
    else
    	fname = cat(1,{SwE.xY.VY.fname}');
    end
    spm_DesRep('DesMtx',SwE.xX,fname,SwE.xsDes)
    fprintf('%30s\n','...done')                                         %-#
end

cd(original_dir); % Change back dir
fprintf('Done\n')