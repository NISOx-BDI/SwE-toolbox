function varargout = swe_list(varargin)
% Display an analysis of an SwE parametric map
% =========================================================================
% FORMAT TabDat = swe_list('List',xSwE,hReg,[Num,Dis,Str])
%    Summary list of local maxima for entire volume of interest
%    xSwE - Structure containing data (format as below)
%    hReg   - Handle of caller (not used)
%    Num    - number of maxima per cluster [3]
%    Dis    - distance among clusters {mm} [8]
%    Str    - header string
%
% FORMAT TabDat = swe_list('ListCluster',xSwE,hReg,[Num,Dis,Str])
%    List of local maxima for a single suprathreshold cluster
%    xSwE - Structure containing data (format as below)
%    hReg   - Handle of caller (not used)
%    Num    - number of maxima per cluster [3]
%    Dis    - distance among clusters {mm} [8]
%    Str    - header string
%
% FORMAT swe_list('TxtList',TabDat,c)
%    Prints a tab-delimited text version of the table
%    TabDat - Structure containing table data (format as below)
%    c      - Column of table data to start text table at
%          (E.g. c=3 doesn't print set-level results contained in columns 1
%           & 2)
%
% FORMAT swe_list('SetCoords',xyz,hAx,hReg)
%    Highlighting of table co-ordinates (used by results section registry)
%    xyz    - 3-vector of new co-ordinate
%    hAx    - table axis (the registry object for tables)
%    hReg   - Handle of caller (not used)
% -------------------------------------------------------------------------
% xSwE   - structure containing SPM, distribution & filtering details
%        - required fields are:
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests
% .STAT  - distribution {Z, T, X or F}
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {voxels}
% .XYZ   - location of voxels {voxel coords}
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}
% .M     - voxels - > mm matrix
% .VOX   - voxel dimensions {mm}
% .DIM   - image dimensions {voxels}
% .units - space units
% .VRpv  - filehandle - Resels per voxel
% .Ps    - uncorrected P values in searched volume (for voxel FDR)
% .Pp    - uncorrected P values of peaks (for peak FDR)
% .Pc    - uncorrected P values of cluster extents (for cluster FDR)
% .uc    - 0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
%          Note: .Ps, .Pp, .Pc and .uc may be specified as empty, i.e. []
% .thresDesc - description of height threshold (string)
%
% (see spm_getSPM.m for further details of xSwE structures)
%
% hReg   - Handle of results section XYZ registry (see swe_results_ui.m)
%
% Num    - number of maxima per cluster [3]
% Dis    - distance among clusters {mm} [8]
% Str    - header string
%
% TabDat - Structure containing table data
%        - fields are
% .tit   - table title (string)
% .hdr   - table header (2x11 cell array)
% .fmt   - fprintf format strings for table data (1x11 cell array)
% .str   - table filtering note (string)
% .ftr   - table footnote information (5x2 cell array)
% .dat   - table data (Nx11 cell array)
%
%                           ----------------
%
% FORMAT swe_list('TxtList',TabDat,c)
% Prints a tab-delimited text version of the table
% TabDat - Structure containing table data (format as above)
% c      - Column of table data to start text table at
%          (E.g. c=3 doesn't print set-level results contained in columns 1 & 2)
%                           ----------------
%
% FORMAT swe_list('SetCoords',xyz,hAx,hReg)
% Highlighting of table co-ordinates (used by results section registry)
% xyz    - 3-vector of new co-ordinate
% hAx    - table axis (the registry object for tables)
% hReg   - Handle of caller (not used)
%__________________________________________________________________________
%
% swe_list characterizes SwE data (thresholded at u and k) in terms of
% excursion sets (a collection of face, edge and vertex connected subsets
% or clusters).  The corrected significance of the results are based on
% voxel-level inferences assuming Gaussian errors, uncorrected or 
% FDR-corrected, or, with wild bootstrap, voxel-level uncorrected, FDR-
% corrected and FWE-corrected, and cluster-level FWE-corrected.
%
% The p values are based on the probability of obtaining c, or more,
% clusters of k, or more, resels above u, in the volume S analysed =
% P(u,k,c).  For specified thresholds u, k, the set-level inference is
% based on the observed number of clusters C, = P(u,k,C).  For each
% cluster of size K the cluster-level inference is based on P(u,K,1)
% and for each voxel (or selected maxima) of height U, in that cluster,
% the voxel-level inference is based on P(U,0,1).  All three levels of
% inference are supported with a tabular presentation of the p values
% and the underlying statistic:
%
% Set-level     - c    = number of suprathreshold clusters
%               - P    = prob(c or more clusters in the search volume)
%
% Cluster-level - k    = number of voxels in this cluster
%               - Pc   = prob(k or more voxels in the search volume)
%               - Pu   = prob(k or more voxels in a cluster)
%               - Qc   = lowest FDR bound for which this cluster would be
%                        declared positive
%
% Peak-level    - T/F  = Statistic upon which the SwE data is based
%               - Ze   = The equivalent Z score - prob(Z > Ze) = prob(t > T)
%               - Pc   = prob(Ze or higher in the search volume)
%               - Qp   = lowest FDR bound for which this peak would be
%                        declared positive
%               - Pu   = prob(Ze or higher at that voxel)
%
% Voxel-level   - Qu   = Expd(Prop of false positives among voxels >= Ze)
%
% x,y,z (mm)    - Coordinates of the voxel
%
% Note: For wild bootstrap settings the following will not be available:
%       - Cluster-level: Pu and Qc
%       - Set-level: P
%
% The table is grouped by regions and sorted on the Ze-variate of the
% primary maxima.  Ze-variates (based on the uncorrected p value) are the
% Z score equivalent of the statistic. Volumes are expressed in voxels.
%
% Clicking on values in the table returns the value to the MATLAB
% workspace. In addition, clicking on the co-ordinates jumps the
% results section cursor to that location. The table has a context menu
% (obtained by right-clicking in the background of the table),
% providing options to print the current table as a text table, or to
% extract the table data to the MATLAB workspace.
%
%_________________________________________________________________________
% Version Info:  $Format:%ci$ $Format:%h$


%==========================================================================
switch lower(varargin{1}), case 'list'                               %-List
%==========================================================================
% FORMAT TabDat = swe_list('List',xSwE,hReg,[Num,Dis,Str])

    %-Parse arguments
    %----------------------------------------------------------------------
    if nargin < 2, error('Not enough input arguments.'); end
    if nargin < 3, hReg = []; else  hReg = varargin{3};  end
    xSwE = varargin{2};

    %-Extract results table and display it
    %----------------------------------------------------------------------
    spm('Pointer','Watch')

    TabDat = swe_list('Table',xSwE,varargin{4:end});

    swe_list('Display',TabDat,hReg);

    spm('Pointer','Arrow')

    %-Return TabDat structure
    %----------------------------------------------------------------------
    varargout = { TabDat };


%==========================================================================
case 'table'                                                        %-Table
%==========================================================================
    % FORMAT TabDat = swe_list('table',xSwE,[Num,Dis,Str])

    %-Parse arguments
    %----------------------------------------------------------------------
    if nargin < 2, error('Not enough input arguments.'); end
    xSwE = varargin{2};

    %-Get number of maxima per cluster to be reported
    %----------------------------------------------------------------------
    if length(varargin) > 2, Num = varargin{3}; else Num = spm_get_defaults('stats.results.volume.nbmax'); end

    %-Get minimum distance among clusters (mm) to be reported
    %----------------------------------------------------------------------
    if length(varargin) > 3, Dis = varargin{4}; else Dis = spm_get_defaults('stats.results.volume.distmin'); end

    %-Get header string
    %----------------------------------------------------------------------
    if length(varargin) > 4 && ~isempty(varargin{5})
        Title = varargin{5};
    else
        if xSwE.STAT ~= 'P'
            Title = 'p-values adjusted for search volume';
        else
            Title = 'Posterior Probabilities';
        end
    end

    %-Extract data from xSwE
    %----------------------------------------------------------------------
    isCifti   = xSwE.isCifti;
    S         = xSwE.S;
    VOX       = xSwE.VOX;
    DIM       = xSwE.DIM;
    M         = xSwE.M;
    XYZ       = xSwE.XYZ;
    Z         = xSwE.Z;
%     VRpv      = xSwE.VRpv;
    n         = xSwE.n;
    STAT      = xSwE.STAT;
    switch STAT
        case 'T'
            STATe = 'Z';
        case 'F'
            STATe = 'X';
    end
    u         = xSwE.u;
    k         = xSwE.k;
    try, uc   = xSwE.uc; end
    try, QPs  = xSwE.Ps; end
    try, QPp  = xSwE.Pp; end
    try, QPc  = xSwE.Pc; end

    % For WB analyses we have already calculated the information for the
    % table and footer. We just need to read it in.
    if xSwE.WB
        VspmUncP = swe_data_read(xSwE.VspmUncP);
        VspmFDRP = swe_data_read(xSwE.VspmFDRP);
        VspmFWEP = swe_data_read(xSwE.VspmFWEP);
        % If the user didn't originally select clusterwise inference,
        % clusterwise FWEP values will not have been calculated.
        if isfield(xSwE, 'VspmFWEP_clus')
            VspmFWEP_clus = swe_data_read(xSwE.VspmFWEP_clus);
        else
            VspmFWEP_clus = [];
        end
        if isfield(xSwE, 'VspmFWEP_clusnorm')
            VspmFWEP_clusnorm = swe_data_read(xSwE.VspmFWEP_clusnorm);
        else
            VspmFWEP_clusnorm = [];
        end
    end

%     if STAT~='P'
%         R     = full(xSwE.R);
%         FWHM  = full(xSwE.FWHM);
%     end
    try
        units = xSwE.units;
    catch
        units = {'mm' 'mm' 'mm'};
    end
    units{1}  = [units{1} ' '];
    units{2}  = [units{2} ' '];

    if ~spm_mesh_detect(xSwE.Vspm)
        DIM   = DIM > 1;                  % non-empty dimensions
        strDataType = 'voxels';
    else
        DIM   = true(1,3);
        strDataType = 'vertices';
    end
    D         = sum(DIM);             % highest dimension
    VOX       = VOX(DIM);             % scaling

%     if STAT ~= 'P'
%         FWHM  = FWHM(DIM);            % Full width at half max
%         FWmm  = FWHM.*VOX;            % FWHM {units}
%         V2R   = 1/prod(FWHM);         % voxels to resels
%         k     = k*V2R;                % extent threshold in resels
%         R     = R(1:(D + 1));         % eliminate null resel counts
%     end
        try, QPs = sort(QPs(:)); end  % Needed for voxel   FDR
        try, QPp = sort(QPp(:)); end  % Needed for peak    FDR
        try, QPc = sort(QPc(:)); end  % Needed for cluster FDR
    % Choose between voxel-wise and topological FDR
    %----------------------------------------------------------------------
    topoFDR = false; %to be checked

    %-Tolerance for p-value underflow, when computing equivalent Z's
    %----------------------------------------------------------------------
    tol = eps*10;

    %-Table Headers
    %----------------------------------------------------------------------
    TabDat.tit = Title;

    % If we are doing a clusterwise/voxelwise analysis the header is the
    % normal SPM header.
    if isCifti
      additionalField = {'brain structure','label','label'};
      nColTable = 12;
    else
      additionalField = {};
      nColTable = 11;
    end

    if ~isfield(xSwE, 'TFCEanaly') || ~xSwE.TFCEanaly
        TabDat.hdr = {...
            'set',      'p',            '\itp';...
            'set',      'c',            '\itc';...
            'cluster',  'p(FWE-corr)',  '\itp\rm_{FWE-corr}';...
            % 'cluster',  'p(FDR-corr)',  '\itq\rm_{FDR-corr}';...
            'cluster',  'equivk',       '\itk\rm_E';...
            'cluster',  'equivkArea',       '\it{area}';...
            'cluster',  'equivkZ',       '\itk\rm_{Z}';...
            % 'cluster',  'p(unc)',       '\itp\rm_{uncorr}';...
            'peak',     'p(FWE-corr)',  '\itp\rm_{FWE-corr}';...
            'peak',     'p(FDR-corr)',  '\itq\rm_{FDR-corr}';...
            'peak',      STATe,         sprintf('\\it%c',STATe);...
            'peak',     'p(unc)',       '\itp\rm_{uncorr}';...
            ' ',         'x,y,z {mm}',   [units{:}];...
            additionalField{:} }';
    % Otherwise we need a TFCE section in the table instead of a cluster
    % level section.
    else
        TabDat.hdr = {...
            'set',      'p',            '\itp';...
            'set',      'c',            '\itc';...
            'TFCE',     'p(FWE-corr)',  '\itp\rm_{FWE-corr}';...
            'TFCE',     '',             '';...
            'TFCE',     'equivk',       '\itk\rm_E';...
            'TFCE',     '',             '';...
            'peak',     'p(FWE-corr)',  '\itp\rm_{FWE-corr}';...
            'peak',     'p(FDR-corr)',  '\itq\rm_{FDR-corr}';...
            'peak',      STATe,         sprintf('\\it%c',STATe);...
            'peak',     'p(unc)',       '\itp\rm_{uncorr}';...
            ' ',         'x,y,z {mm}',   [units{:}];...
            additionalField{:} }';
    end

    %-Coordinate Precisions
    %----------------------------------------------------------------------
    if isempty(xSwE.XYZmm) || isCifti % empty results or cifti
        xyzfmt = '%3.0f %3.0f %3.0f';
        voxfmt = '%1.1f %1.1f %1.1f';
    elseif ~any(strcmp(units{3},{'mm',''})) % 2D data
        xyzfmt = '%3.0f %3.0f %3.0f';
        voxfmt = '%1.1f %1.1f %1.1f';
    else % 3D data, work out best precision based on voxel sizes and FOV
        xyzsgn = min(xSwE.XYZmm(DIM,:),[],2) < 0;
        xyzexp = max(floor(log10(max(abs(xSwE.XYZmm(DIM,:)),[],2)))+(max(abs(xSwE.XYZmm(DIM,:)),[],2) >= 1),0);
        voxexp = floor(log10(abs(VOX')))+(abs(VOX') >= 1);
        xyzdec = max(-voxexp,0);
        voxdec = max(-voxexp,1);
        xyzwdt = xyzsgn+xyzexp+(xyzdec>0)+xyzdec;
        voxwdt = max(voxexp,0)+(voxdec>0)+voxdec;
        tmpfmt = cell(size(xyzwdt));
        for i = 1:numel(xyzwdt)
            tmpfmt{i} = sprintf('%%%d.%df ', xyzwdt(i), xyzdec(i));
        end
        xyzfmt = [tmpfmt{:}];
        tmpfmt = cell(size(voxwdt));
        for i = 1:numel(voxwdt)
            tmpfmt{i} = sprintf('%%%d.%df ', voxwdt(i), voxdec(i));
        end
        voxfmt = [tmpfmt{:}];
    end
    TabDat.fmt = {'%-0.3f','%g',...                            %-Set
        '%0.3f', '%0.0f','%0.2f','%0.3f',...                   %-Cluster
        '%0.3f', '%0.3f', '%6.2f', '%0.3f',...                 %-Peak
        xyzfmt, '%s'};                                         %-XYZ

    %-Table filtering note
    %----------------------------------------------------------------------
    if isinf(Num)
        TabDat.str = sprintf('table shows all local maxima > %.1fmm apart',Dis);
    else
        TabDat.str = sprintf(['table shows %d local maxima ',...
            'more than %.1fmm apart'],Num,Dis);
    end

    %-Footnote with SPM parameters
    %----------------------------------------------------------------------
     if strcmp(STAT, 'T')
         Pz              = 1-spm_Ncdf(u);
         eSTAT = 'Z';
     else
         Pz              = 1-spm_Xcdf(u, 1);
         eSTAT = 'X';
     end

     % Create footer for display.
     TabDat.ftr      = cell(6,2);

     % Number of `extra` lines inserted that don't have to be present in
     % every display
     exlns           = 0;

     % detect whether the WB was done based on Z/X or on T/F using the version number (2.1.1 was the last using T/F)
     if swe_compareVersions(swe('ver'), '2.1.1', '>')
      displaySTAT = eSTAT;
     else
      displaySTAT = STAT;
     end

     if xSwE.WB
        % Record thresholds.
        % For voxel-wise FDR and unc. thresholding, we cannot display a score threshold as it varies per voxel
        td = regexp(xSwE.thresDesc,'p\D+?(?<u>[\.\d]+) \((?<thresDesc>\S+)\)','names');
        if xSwE.infType == 0 % voxel-wise
            if strcmp(td.thresDesc, 'FWE')
              TabDat.ftr{1,1} = ['Threshold: Height ' displaySTAT ' > %0.2f, p <= %0.3f (FWE); Extent k >= %0.0f ' strDataType '.'];
              TabDat.ftr{1,2} = [u, str2num(td.u), k];
            elseif strcmp(td.thresDesc, 'FDR')
              TabDat.ftr{1,1} = ['Threshold: p <= %0.3f (FDR); Extent k >= %0.0f ' strDataType '.'];
              TabDat.ftr{1,2} = [str2num(td.u), k];
            elseif strcmp(td.thresDesc, 'unc.')
              TabDat.ftr{1,1} = ['Threshold: p <= %0.3f (unc.); Extent k >= %0.0f ' strDataType '.'];
              TabDat.ftr{1,2} = [str2num(td.u), k];
            else
              error('Unknown inference type detected')
            end
        elseif xSwE.infType == 1 % cluster-wise
            if strcmp(xSwE.clustWise, 'FWE')
              if strcmpi(xSwE.clusterSizeType, 'Box-Cox norm. k_{Z}')
                TabDat.ftr{1,1} = ['Threshold: Height ' eSTAT ' > %0.2f, p < %0.3f (unc.); k_{Z} > %0.3f, p <= %0.3f (FWE).'];
                TabDat.ftr{1,2} = [u, str2num(td.u), k, xSwE.fwep_c];
              elseif strcmpi(xSwE.clusterSizeType, 'classic k_E')
                TabDat.ftr{1,1} = ['Threshold: Height ' eSTAT ' > %0.2f, p < %0.3f (unc.); Extent k > %0.0f ' strDataType ', p <= %0.3f (FWE).'];
                TabDat.ftr{1,2} = [u, str2num(td.u), k, xSwE.fwep_c];
              else
                error('Unknow type of cluster statistics');
              end
            elseif strcmp(xSwE.clustWise, 'Uncorr')
              TabDat.ftr{1,1} = ['Threshold: p <= %0.3f (unc.); Extent k >= %0.0f ' strDataType '.'];
              TabDat.ftr{1,2} = [str2num(td.u), k];
            else
              error('Unknown inference type detected')
            end
        elseif xSwE.infType == 2 % TFCE
            TabDat.ftr{1,1} = 'Threshold: TFCE %s';
            TabDat.ftr{1,2} = xSwE.thresDesc;
        else
            error('Unknown inference type detected')
        end
        % We need the P uncorrected P values to be in the correct form to
        % use spm_uc_FDR.
        % Make sure to load only the in-mask data to avoid loading zero values
        Ts = swe_data_read(xSwE.VspmUncP, 'xyz', xSwE.XYZ_inMask);
        Ts(isnan(Ts)) = [];
        Ts = 10.^-Ts;
        Ts = sort(Ts(:));

        % Obtain the FDR p 0.05 value.
        FDRp_05 = swe_uc_FDR(0.05,Inf,'P',n,Ts);
        clear Ts

        % Record FWE/FDR/clus FWE p values. (No clus FWE for voxelwise and
        % TFCE analyses)
        if xSwE.infType == 1 && strcmp(xSwE.clustWise, 'FWE')
            TabDat.ftr{2,1} = ...
                 ['vox ' displaySTAT '(5%% FWE): %0.3f, vox P(5%% FDR): %0.3f, clus k(5%% FWE): %0.0f '];
            TabDat.ftr{2,2} = [xSwE.Pfv, FDRp_05, xSwE.Pfc];
        else
            TabDat.ftr{2,1} = ...
                 ['vox ' displaySTAT '(5%% FWE): %0.3f, vox P(5%% FDR): %0.3f'];
            TabDat.ftr{2,2} = [xSwE.Pfv, FDRp_05];
        end
     else
        % Record height thresholds.
        TabDat.ftr{1,1} = ...
        ['Threshold: Height ' eSTAT ' = %0.2f, p = %0.3f; Extent k = %0.0f ' strDataType '.'];
        TabDat.ftr{1,2} = [u,Pz,k];
        % Record FDR p value.
        TabDat.ftr{2,1} = ...
             'vox P(5%% FDR): %0.3f';
        TabDat.ftr{2,2} = swe_uc_FDR(0.05,Inf,'P',n,sort(xSwE.Ps)');

     end

     if xSwE.infType == 1 && strcmp(xSwE.clustWise, 'FWE') && isfield(xSwE, 'boxcoxInfo')
        TabDat.ftr{(3+exlns),1} = 'k_{Z} = 0.6745 (k_{\\lambda} - Q2(k_{\\lambda}^{H0})) / (Q3(k_{\\lambda}^{H0})-Q2(k_{\\lambda}^{H0}))';
        % TabDat.ftr{(3+exlns),1} = 'Null cluster sizes in surfaces: \\lambda_S=%0.2f , \\lambda_V =%0.2f';
        % TabDat.ftr{(3+exlns),1} = 'Box-Cox(Surf): \lambda=%0.2f, mean=%0.2f, std=%0.2f, median=%0.2f, 2(Q3-Q2)=%0.2f';
        TabDat.ftr{(3+exlns),2} = [];
        if isfield(xSwE.boxcoxInfo, 'surfaces') && isfield(xSwE.boxcoxInfo, 'volume')
            tmpSurf = xSwE.boxcoxInfo.surfaces;
            tmpVol = xSwE.boxcoxInfo.volume;
            TabDat.ftr{(4+exlns),1} = 'Box-Cox parameters for cluster sizes under H0: \\lambda(Surfaces)=%0.2f, \\lambda(Volume)=%0.2f';
            TabDat.ftr{(4+exlns),2} = [tmpSurf.lambda, tmpVol.lambda];
        elseif isfield(xSwE.boxcoxInfo, 'surfaces')
            tmpSurf = xSwE.boxcoxInfo.surfaces;
            TabDat.ftr{(4+exlns),1} = 'Box-Cox parameters for cluster sizes under H0: \\lambda(Surfaces)=%0.2f';
            TabDat.ftr{(4+exlns),2} = [tmpSurf.lambda];
        elseif isfield(xSwE.boxcoxInfo, 'volume')
            tmpVol = xSwE.boxcoxInfo.volume;
            TabDat.ftr{(4+exlns),1} = 'Box-Cox parameters for cluster sizes under H0: \\lambda(Volume)=%0.2f';
            TabDat.ftr{(4+exlns),2} = [tmpVol.lambda];
        else
            error('Unknown Box-Cox Info!')
        end
        exlns = exlns + 2;
     end

     % If we have groups display group details in ftr.
     if isfield(xSwE, 'nSubj_g')
         % Record number of subjects per group.
         nSubjString = 'Number of subjects: ';
         for i = 1:length(xSwE.nSubj_g)
             nSubjString = [nSubjString '%0.0f'];
             if i ~= length(xSwE.nSubj_g)
                 nSubjString = [nSubjString ', '];
             else
                 nSubjString = [nSubjString '; '];
             end
         end

         % Record visits per group.
         nVisitsString = 'Number of visits ([Min Max]): ';
         nVisitsNumbers = [];
         for i = 1:length(xSwE.max_nVis_g)
             nVisitsString = [nVisitsString '[%0.0f %0.0f]'];
             if i ~= length(xSwE.max_nVis_g)
                nVisitsString = [nVisitsString ', '];
             end
             nVisitsNumbers = [nVisitsNumbers xSwE.min_nVis_g(i) xSwE.max_nVis_g(i)];
         end
         TabDat.ftr{3+exlns,1} = [nSubjString nVisitsString];
         TabDat.ftr{3+exlns,2} = [xSwE.nSubj_g nVisitsNumbers];

         exlns = exlns + 1;
     end

     % Retrieve edf data
     if isfield(xSwE, 'Vedf')
        edf = swe_data_read(xSwE.Vedf, 'xyz', xSwE.XYZ_inMask);
     else
        edf = xSwE.edf;
     end
     edf(isnan(edf)) = [];

     edf_max = max(edf);
     edf_min = min(edf);
     edf_med = median(edf);

     % Work out range of dof values
     diff = abs(edf_max - edf_min);

     % Work out dofType
     switch xSwE.dofType

         case 0

             dofTypeStr = 'Naive';

         case 1

             dofTypeStr = 'Approx I';

         case 2

             dofTypeStr = 'Approx II';

         case 3

             dofTypeStr = 'Approx III';

         otherwise

             error('Unknown degrees of freedom.')

     end

     % Recording effective Degrees of freedom
     if xSwE.dofType~=0 && diff > 10^-10
        TabDat.ftr{(3+exlns),1}=['Error DF: (' dofTypeStr '): (min) %0.1f, (median) %0.1f, (max) %0.1f'];
        TabDat.ftr{(3+exlns),2}=[edf_min, edf_med, edf_max];
     else
        TabDat.ftr{(3+exlns),1}=['Error DF: (' dofTypeStr '): %0.1f'];
        TabDat.ftr{(3+exlns),2}=edf_med;
     end

    % Record small sample adjustments.
    TabDat.ftr{(4+exlns),1}='Resid. Adj.: %s';
    switch xSwE.SS
        case {0, 1, 2, 3}
            TabDat.ftr{(4+exlns),2} = ['Type ' num2str(xSwE.SS)];
        case {4, 5}
            TabDat.ftr{(4+exlns),2} = ['Type C' num2str(xSwE.SS - 2)];
        otherwise
            error('Unknown SS type')
    end
     % Record contrast degrees of freedom.
     TabDat.ftr{(5+exlns),1} = 'Contrast DF: %0.0f; Number of predictors: %0.0f';
     TabDat.ftr{(5+exlns),2} = [xSwE.df_Con xSwE.nPredict];

     % Record volume.
     if isCifti
        TabDat.ftr{(6+exlns),1} = ...
        ['Surface(s): %0.0f vertices; Volume: %0.0f voxels'];
        TabDat.ftr{(6+exlns),2} = [xSwE.S_surf, xSwE.S_vol];
     elseif spm_mesh_detect(xSwE.Vspm)
        TabDat.ftr{(6+exlns),1} = ...
            ['Surface: %0.0f ' strDataType ''];
        TabDat.ftr{(6+exlns),2} = [S];
     else
        TabDat.ftr{(6+exlns),1} = ...
            ['Volume: %0.0f ' units{:} ' = %0.0f ' strDataType ''];
        TabDat.ftr{(6+exlns),2} = [S*prod(VOX),S];
     end

     % Record voxel sizes.
     if isCifti && numel(xSwE.cifti.volume) > 0
        TabDat.ftr{(7+exlns),1} = ...
            ['Voxel size: %1.1f %1.1f %1.1f mm mm mm'];
        TabDat.ftr{(7+exlns),2} = sqrt(diag(xSwE.cifti.volume.M(1:3,1:3)'*xSwE.cifti.volume.M(1:3,1:3)))';
     elseif isCifti && numel(xSwE.cifti.volume) == 0
        exlns = exlns - 1;
     elseif ~spm_mesh_detect(xSwE.Vspm)
        TabDat.ftr{(7+exlns),1} = ...
            ['Voxel size: ' voxfmt units{:}];
        TabDat.ftr{(7+exlns),2} = VOX;
     else
        exlns = exlns - 1;
     end

     if isfield(xSwE, 'TFCEanaly') && xSwE.TFCEanaly
         TabDat.ftr{(8+exlns),1} = 'TFCE: E=%0.1f, H=%0.1f';
         TabDat.ftr{(8+exlns),2} = [xSwE.TFCE.E, xSwE.TFCE.H];
         exlns = exlns + 1;
     end
     if xSwE.WB

        % Recording number of bootstraps.
        TabDat.ftr{(8+exlns),1}='Bootstrap samples = %0.0f';
        TabDat.ftr{(8+exlns),2}= xSwE.nB;

        exlns = exlns + 1;

    end
    %-Characterize excursion set in terms of maxima
    % (sorted on Z values and grouped by regions)
    %----------------------------------------------------------------------
    if isempty(Z)
        TabDat.dat = cell(0,nColTable);
        varargout  = {TabDat};
        return
    end

    %-Workaround in spm_max for conjunctions with negative thresholds
    %----------------------------------------------------------------------
    minz           = abs(min(min(Z)));
    Z              = 1 + minz + Z;

    if xSwE.infType == 1 && strcmp(xSwE.clustWise, 'FWE') && isfield(xSwE, 'boxcoxInfo')
        boxcoxInfo = xSwE.boxcoxInfo;
    else
        boxcoxInfo = [];
    end

    if isCifti
        [N, N_area, N_boxcox, Z, XYZ, A, L, XYZmm, brainStructureShortLabels] = ...
            swe_cifti_max(Z,XYZ(1,:), xSwE.cifti, boxcoxInfo);
    elseif ~spm_mesh_detect(xSwE.Vspm)
        [N, N_boxcox, Z, XYZ, A, L]  = swe_max(Z, XYZ, boxcoxInfo);
        N_area = [];
    else
        [N, N_area, N_boxcox, Z, XYZ, A, L]  = swe_mesh_max(Z, XYZ, xSwE.G, boxcoxInfo);
    end
    Z              = Z - minz - 1;

    %-Convert cluster sizes from voxels (N) to resels (K)
    %----------------------------------------------------------------------
    c              = max(A);                           %-Number of clusters
    NONSTAT        = spm_get_defaults('stats.rft.nonstat');

    %-Convert maxima locations from voxels to mm
    %----------------------------------------------------------------------
    if isCifti
      % nothing as it was done in swe_cifti_max above
    elseif spm_mesh_detect(xSwE.Vspm)
        XYZmm = xSwE.G.vertices(XYZ(1,:),:)';
    else
        XYZmm = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];
    end

    %-Set-level p values {c} - do not display if reporting a single cluster
    %----------------------------------------------------------------------
%     if STAT ~= 'P'
%         Pc     = spm_P(c,k,u,df,STAT,R,n,S);            %-Set-level p-value
%     else
        Pc     = [];
%     end

    TabDat.dat = {Pc,c};
    TabLin     = 1;

    %-Cluster and local maxima p-values & statistics
    %----------------------------------------------------------------------
    while numel(find(isfinite(Z)))

        %-Find largest remaining local maximum
        %------------------------------------------------------------------
        [U,i]  = max(Z);            %-largest maxima
        mj     = find(A == A(i));   %-maxima in cluster

        if ~xSwE.WB
            switch STATe
                case 'Z'
                try
                    Pz      = normcdf(-U);
                catch
                    Pz      = spm_Ncdf(-U);
                end
                case 'X'
                try
                    Pz      = 1-chi2cdf(U,1);
                catch
                    Pz      = 1-spm_Xcdf(U,1);
                end
            end
        else
            Pz = 10.^-VspmUncP(XYZ(1,i),XYZ(2,i),XYZ(3,i));
        end

        % If we are not running a wild bootstrap or we are doing a
        % small volume correction we need to calculate the FDR P value
        % and leave the other values blank.
        if ~xSwE.WB || isfield(xSwE,'svc')
            Pu      = [];
            Pk      = [];
            Pn      = [];
            Qc      = [];
            Qp      = [];
            Qu      = spm_P_FDR(U,[1 1],STATe,n,QPs);    % voxel FDR-corrected
            ws      = warning('off','SPM:outOfRangeNormal');
            warning(ws);
        % If we are running a wild bootstrap we only need to read in
        % results we calculated earlier.
        else
            Pu      = 10.^-VspmFWEP(XYZ(1,i),XYZ(2,i),XYZ(3,i));
            Qu      = 10.^-VspmFDRP(XYZ(1,i),XYZ(2,i),XYZ(3,i));
            Pn      = [];
            Qc      = [];
            Qp      = [];
            ws      = warning('off','SPM:outOfRangeNormal');
            warning(ws);

            if xSwE.infType == 1 && strcmp(xSwE.clustWise, 'FWE') % only for FWER clusterwise WB
                if strcmpi(xSwE.clusterSizeType, 'Box-Cox norm. k_{Z}')
                    Pk  = 10.^-VspmFWEP_clusnorm(XYZ(1,i),XYZ(2,i),XYZ(3,i));
                elseif strcmpi(xSwE.clusterSizeType, 'classic k_E')
                    Pk  = 10.^-VspmFWEP_clus(XYZ(1,i),XYZ(2,i),XYZ(3,i));
                else
                    error('Unknow type of cluster statistics');
                end
                % It is possible to get the results window to display
                % details about voxels that were thresholded out by the
                % cluster-forming threshold used for the wild bootstrap.
                % These regions will have NaN for the cluster FWE P-value
                % when they should have one. So the below is necessary:
                Pk(isnan(Pk)) = 1;

            elseif xSwE.TFCEanaly

                % Get coordinates of all voxels in the current cluster.
                currentClus = find(A == A(i));
                XYZ_clus = XYZ(:, currentClus);

                % Read in all TFCE FWE P values in this cluser
                tfp = 10.^-swe_data_read(xSwE.VspmTFCEFWEP, 'xyz', XYZ_clus);

                % Record the minimum TFCE FWE P value in said cluster.
                Pk  = min(tfp);
            else
                Pk  = [];
            end

        end

        if i > numel(N_area) % means that this is for volume or there is no area info
            N_area_tmp = [];
        else
            N_area_tmp = N_area(i);
        end
        if i > numel(N_boxcox)
            N_boxcox_tmp = [];
        else
            N_boxcox_tmp = N_boxcox(i);
        end
        [TabDat.dat{TabLin,3:11}] = deal(Pk, N(i), N_area_tmp, N_boxcox_tmp,Pu,Qu,U,Pz,XYZmm(:,i));

        if isCifti
          [TabDat.dat{TabLin, 12}] = char(brainStructureShortLabels(i));
        end
        TabLin = TabLin + 1;

        %-Print Num secondary maxima (> Dis mm apart)
        %------------------------------------------------------------------
        [l,q] = sort(-Z(mj));                              % sort on Z value
        D     = i;
        for i = 1:length(q)
            d = mj(q(i));
            if min(sqrt(sum((XYZmm(:,D)-repmat(XYZmm(:,d),1,size(D,2))).^2)))>Dis
                if length(D) < Num
                    % voxel-level p values {Z}
                    %------------------------------------------------------
%                     if STAT ~= 'P'
%                         Pz     = spm_P(1,0,Z(d),df,STAT,1,n,S);
%                         Pu     = spm_P(1,0,Z(d),df,STAT,R,n,S);
%                         if topoFDR
%                             Qp = spm_P_peakFDR(Z(d),df,STAT,R,n,u,QPp);
%                             Qu = [];
%                         else
%                             Qu = spm_P_FDR(Z(d),df,STATe,n,QPs);
%                             Qp = [];
%                         end
%                         if Pz < tol
%                             Ze = Inf;
%                         else
%                             Ze = swe_invNcdf(1 - Pz);
%                         end
%                     else

                    % If we are not running a wild bootstrap or if we are
                    % doing a small volume correction we need to calculate
                    % the FDR P value and leave the other values blank.
                    if ~xSwE.WB || isfield(xSwE,'svc')
                        Pz     = spm_Ncdf(-Z(d));
                        Pu     = [];
                        Qu     = [];
                        Qp     = [];
                        Qu      = spm_P_FDR(Z(d),[1 1],STATe,n,QPs);    % voxel FDR-corrected
                        ws     = warning('off','SPM:outOfRangeNormal');
                        Ze     = swe_invNcdf(Z(d));
                        warning(ws);
%                     end

                    % If we are running a wild bootstrap we only need to read in
                    % results we calculated earlier.
                    else

                        Pz      = 10.^-VspmUncP(XYZ(1,d),XYZ(2,d),XYZ(3,d));
                        Pu      = 10.^-VspmFWEP(XYZ(1,d),XYZ(2,d),XYZ(3,d));
                        Qu      = 10.^-VspmFDRP(XYZ(1,d),XYZ(2,d),XYZ(3,d));
                        ws     = warning('off','SPM:outOfRangeNormal');
                        Ze     = swe_invNcdf(Z(d));
                        warning(ws);

                    end

                    D     = [D d];
                    if topoFDR
                    [TabDat.dat{TabLin,7:11}] = ...
                        deal(Pu,Qp,Z(d),Pz,XYZmm(:,d));
                    else
                    [TabDat.dat{TabLin,7:11}] = ...
                        deal(Pu,Qu,Z(d),Pz,XYZmm(:,d));
                    end
                    TabLin = TabLin+1;
                end
            end
        end
        Z(mj) = NaN;     % Set local maxima to NaN
    end

    varargout = {TabDat};

    %======================================================================
    case 'display'                       %-Display table in Graphics window
    %======================================================================
    % FORMAT swe_list('display',TabDat,hReg)

    %-Parse arguments
    %----------------------------------------------------------------------
    if nargin < 2, error('Not enough input arguments.');
    else           TabDat = varargin{2}; end
    if nargin < 3, hReg = []; else hReg = varargin{3}; end

    isCifti = (size(TabDat.hdr,2) == 12);
    if isCifti
      scalingFactor = 0.9;
      nColTable = 12;
    else
      scalingFactor = 1;
      nColTable = 11;
    end

    %-Get current location (to highlight selected voxel in table)
    %----------------------------------------------------------------------
    xyzmm = swe_results_ui('GetCoords');

    %-Setup Graphics panel
    %----------------------------------------------------------------------
    Fgraph = spm_figure('FindWin','Satellite');
    if ~isempty(Fgraph)
        spm_figure('Focus', Fgraph);
        ht = 0.85; bot = 0.14;
    else
        Fgraph = spm_figure('GetWin','Graphics');
        ht = 0.4; bot = 0.1;
    end
    swe_results_ui('Clear',Fgraph)

    FS     = spm('FontSizes');           %-Scaled font sizes
    PF     = spm_platform('fonts');      %-Font names (for this platform)

    %-Table axes & Title
    %----------------------------------------------------------------------
    hAx   = axes('Parent',Fgraph,...
                 'Position',[0.025 bot 0.9 ht],...
                 'DefaultTextFontSize',FS(8),...
                 'DefaultTextInterpreter','Tex',...
                 'DefaultTextVerticalAlignment','Baseline',...
                 'Tag','SPMList',...
                 'Units','points',...
                 'Visible','off');

    try
        hRotate3d = rotate3d(Fgraph);
        setAllowAxesRotate(hRotate3d, hAx, false);
    end

    AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
    dy    = FS(9);
    y     = floor(AxPos(4)) - dy;

    text(0,y,['Statistics:  \it\fontsize{',num2str(FS(9)),'}',TabDat.tit],...
              'FontSize',FS(11),'FontWeight','Bold');   y = y - dy/2;
    line([0 1],[y y],'LineWidth',3,'Color','r'),        y = y - 9*dy/8;

    %-Display table header
    %----------------------------------------------------------------------
    set(gca,'DefaultTextFontName',PF.helvetica,'DefaultTextFontSize',FS(8))

    Hs = []; Hc = []; Hp = [];
    h  = text(0.01*scalingFactor,y, [TabDat.hdr{1,1} '-level'],'FontSize',FS(9));
    h  = line([0,0.11]*scalingFactor,[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
    h  = text(0.02*scalingFactor,y-9*dy/8,    TabDat.hdr{3,1});              Hs = [Hs,h];
    h  = text(0.08*scalingFactor,y-9*dy/8,    TabDat.hdr{3,2});              Hs = [Hs,h];

    text(0.22*scalingFactor,y, [TabDat.hdr{1,3} '-level'],'FontSize',FS(9));
    line([0.14,0.44]*scalingFactor,[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
    h  = text(0.15*scalingFactor,y-9*dy/8,    TabDat.hdr{3,3});              Hc = [Hc,h];
    h  = text(0.24*scalingFactor,y-9*dy/8,    TabDat.hdr{3,4});              Hc = [Hc,h];
    h  = text(0.31*scalingFactor,y-9*dy/8,    TabDat.hdr{3,5});              Hc = [Hc,h];
    h  = text(0.39*scalingFactor,y-9*dy/8,    TabDat.hdr{3,6});              Hc = [Hc,h];

    text(0.59*scalingFactor,y, [TabDat.hdr{1,8} '-level'],'FontSize',FS(9));
    line([0.48,0.80]*scalingFactor,[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
    h  = text(0.49*scalingFactor,y-9*dy/8,    TabDat.hdr{3,7});              Hp = [Hp,h];
    h  = text(0.58*scalingFactor,y-9*dy/8,    TabDat.hdr{3,8});              Hp = [Hp,h];
    h  = text(0.67*scalingFactor,y-9*dy/8,    TabDat.hdr{3,9});              Hp = [Hp,h];
    h  = text(0.74*scalingFactor,y-9*dy/8,    TabDat.hdr{3,10});             Hp = [Hp,h];

    text(0.845*scalingFactor,y - dy/2,TabDat.hdr{3,11},'Fontsize',FS(8));

    if isCifti
      text(0.88,y - dy/2,TabDat.hdr{1,12},'Fontsize',FS(8));
    end

    %-Move to next vertical position marker
    %----------------------------------------------------------------------
    y     = y - 7*dy/4;
    line([0 1],[y y],'LineWidth',1,'Color','r')
    y     = y - 5*dy/4;
    y0    = y;

    %-Table filtering note
    %----------------------------------------------------------------------
    text(0.5,4,TabDat.str,'HorizontalAlignment','Center',...
        'FontName',PF.helvetica,'FontSize',FS(8),'FontAngle','Italic')

    %-Footnote with SPM parameters (if classical inference)
    %----------------------------------------------------------------------
    line([0 1],[0 0],'LineWidth',1,'Color','r')
    if ~isempty(TabDat.ftr)
        set(gca,'DefaultTextFontName',PF.helvetica,...
            'DefaultTextInterpreter','Tex','DefaultTextFontSize',FS(8))

        fx = repmat([0 0.645],ceil(size(TabDat.ftr,1)/2),1);
        fy = repmat((1:ceil(size(TabDat.ftr,1)/2))',1,2);
        for i=1:size(TabDat.ftr,1)
            text(fx(i),-fy(i)*dy,sprintf(TabDat.ftr{i,1},TabDat.ftr{i,2}),...
                'UserData',TabDat.ftr{i,2},...
                'ButtonDownFcn','get(gcbo,''UserData'')');
        end
    end

    %-Characterize excursion set in terms of maxima
    % (sorted on Z values and grouped by regions)
    %======================================================================
    if isempty(TabDat.dat)
        text(0.5,y-6*dy,'no suprathreshold clusters',...
            'HorizontalAlignment','Center',...
            'FontAngle','Italic','FontWeight','Bold',...
            'FontSize',FS(16),'Color',[1,1,1]*.5);
        return
    end

    %-Table proper
    %======================================================================

    %-Column Locations
    %----------------------------------------------------------------------
    tCol = [ 0.01      0.08 ...                                %-Set
             0.15      0.24      0.31      0.39 ...            %-Cluster
             0.49      0.58      0.65      0.74 ...            %-Peak
             0.84 ...                                          %-XYZ
             0.88/scalingFactor ] * scalingFactor;                                            %-Brain structure labels

    %-Pagination variables
    %----------------------------------------------------------------------
    hPage = [];
    set(gca,'DefaultTextFontName',PF.courier,'DefaultTextFontSize',FS(7));

    %-Set-level p values {c} - do not display if reporting a single cluster
    %----------------------------------------------------------------------
     if isempty(TabDat.dat{1,1}) && isempty(TabDat.dat{1,2}) % Pc
        set(Hs,'Visible','off');
     end

    if TabDat.dat{1,2} >= 1 % c
        h     = text(tCol(1),y,sprintf(TabDat.fmt{1},TabDat.dat{1,1}),...
                    'FontWeight','Bold', 'UserData',TabDat.dat{1,1},...
                    'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(2),y,sprintf(TabDat.fmt{2},TabDat.dat{1,2}),...
                    'FontWeight','Bold', 'UserData',TabDat.dat{1,2},...
                    'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
    else
        set(Hs,'Visible','off');
    end

    %-Cluster and local maxima p-values & statistics
    %----------------------------------------------------------------------
    HlistXYZ = [];
    for i=1:size(TabDat.dat,1)

        %-Paginate if necessary
        %------------------------------------------------------------------
        if y < dy
%             h = text(0.5,-5*dy,...
%                 sprintf('Page %d',spm_figure('#page',Fgraph)),...
%                         'FontName',PF.helvetica,'FontAngle','Italic',...
%                         'FontSize',FS(8));
            spm_figure('NewPage',hPage)
            hPage = [];
            y     = y0;
        end

        %-Print cluster and maximum peak-level p values
        %------------------------------------------------------------------
        if  ~isempty(TabDat.dat{i,4}), fw = 'Bold'; else fw = 'Normal'; end

        for k=3:nColTable
          if k < 11
            h     = text(tCol(k),y,sprintf(TabDat.fmt{k},TabDat.dat{i,k}),...
                        'FontWeight',fw,...
                        'UserData',TabDat.dat{i,k},...
                        'ButtonDownFcn','get(gcbo,''UserData'')', 'Interpreter', 'none');
            hPage = [hPage, h];
          elseif k == 12
            BDFcn  = ['shortLabel = get(gcbo,''UserData'');',...
                      'if strcmp(shortLabel, ''S_L''),',...
                      'longLabels = ''CIFTI_STRUCTURE_CORTEX_LEFT'';',...
                      'elseif strcmp(shortLabel, ''S_R''),',...
                      'longLabels = ''CIFTI_STRUCTURE_CORTEX_RIGHT'';',...
                      'else,',...
                      'longLabels = swe_cifti_convertVolLabels(extractAfter(shortLabel, ''V_''), false);',...
                      'end,',...
                      'sprintf(''%s: %s'', shortLabel, longLabels)'];
            h     = text(tCol(k),y,sprintf(TabDat.fmt{k},TabDat.dat{i,k}),...
                        'FontWeight',fw,...
                        'UserData',TabDat.dat{i,k},...
                        'ButtonDownFcn', BDFcn, 'Interpreter', 'none');
            hPage = [hPage, h];
          end
        end

        % Specifically changed so it properly finds hMIPax
        %------------------------------------------------------------------
        tXYZmm = TabDat.dat{i,11};
        BDFcn  = [...
            'spm_mip_ui(''SetCoords'',get(gcbo,''UserData''),',...
            'findobj(''tag'',''hMIPax''));'];
        BDFcn = 'spm_XYZreg(''SetCoords'',get(gcbo,''UserData''),hReg,1);';
        h      = text(tCol(11),y,sprintf(TabDat.fmt{11},tXYZmm),...
            'FontWeight',fw,...
            'Tag','ListXYZ',...
            'ButtonDownFcn',BDFcn,...
            'Interruptible','off',...
            'BusyAction','Cancel',...
            'UserData',tXYZmm);

        HlistXYZ = [HlistXYZ, h];
        if spm_XYZreg('Edist',xyzmm,tXYZmm)<eps && ~isempty(hReg)
            set(h,'Color','r')
        end
        hPage  = [hPage, h];

        y      = y - dy;
    end

    %-Number and register last page (if paginated)
    %----------------------------------------------------------------------
    if spm_figure('#page',Fgraph)>1
%         h = text(0.5,-5*dy,sprintf('Page %d/%d',spm_figure('#page',Fgraph)*[1,1]),...
%             'FontName',PF.helvetica,'FontSize',FS(8),'FontAngle','Italic');
        spm_figure('NewPage',hPage)
    end

    %-End: Store TabDat in UserData of context menu
    %======================================================================
    h = uicontextmenu('Tag','TabDat','UserData',TabDat);
    set(gca,'UIContextMenu',h,...
        'Visible','on',...
        'XColor','w','YColor','w')
    uimenu(h,'Label','Print text table',...
        'CallBack',...
        'swe_list(''txtlist'',get(get(gcbo,''Parent''),''UserData''),3)',...
        'Interruptible','off','BusyAction','Cancel');
    uimenu(h,'Label','Extract table data structure',...
        'CallBack','TabDat=get(get(gcbo,''Parent''),''UserData'')',...
        'Interruptible','off','BusyAction','Cancel');
    if ispc
        uimenu(h,'Label','Export to Excel',...
        'CallBack',...
        'swe_list(''xlslist'',get(get(gcbo,''Parent''),''UserData''))',...
        'Interruptible','off','BusyAction','Cancel');
    end
    uimenu(h,'Label','Export to CSV file',...
        'CallBack',...
        'swe_list(''csvlist'',get(get(gcbo,''Parent''),''UserData''))',...
        'Interruptible','off','BusyAction','Cancel');

    %-Setup registry
    %----------------------------------------------------------------------
    set(hAx,'UserData',struct('hReg',hReg,'HlistXYZ',HlistXYZ))
    spm_XYZreg('Add2Reg',hReg,hAx,'swe_list');

    varargout = {};

    %======================================================================
    case 'listcluster'                      %-List for current cluster only
    %======================================================================
    % FORMAT TabDat = swe_list('ListCluster',xSwE,hReg,[Num,Dis,Str])

        %-Parse arguments
        %------------------------------------------------------------------
        if nargin < 2, error('Not enough input arguments.'); end
        if nargin < 3, hReg = []; else hReg = varargin{3};   end
        xSwE = varargin{2};

        if isfield(xSwE,'G')
            warning('"current cluster" option not implemented for meshes.');
            varargout = { evalin('base','TabDat') };
            return;
        end

        %-Get number of maxima per cluster to be reported
        %------------------------------------------------------------------
        if nargin < 4, Num = spm_get_defaults('stats.results.svc.nbmax'); else Num = varargin{4}; end

        %-Get minimum distance among clusters (mm) to be reported
        %------------------------------------------------------------------
        if nargin < 5, Dis = spm_get_defaults('stats.results.svc.distmin'); else Dis = varargin{5}; end

        %-Get header string
        %------------------------------------------------------------------
        if nargin < 6, Str = ''; else Str = varargin{6}; end

        %-If there are suprathreshold voxels, filter out all but current cluster
        %------------------------------------------------------------------
        if ~isempty(xSwE.Z)

            %-Jump to voxel nearest current location
            %--------------------------------------------------------------
            [xyzmm,i] = spm_XYZreg('NearestXYZ',...
                swe_results_ui('GetCoords'),xSwE.XYZmm);
            swe_results_ui('SetCoords',xSwE.XYZmm(:,i));

            %-Find selected cluster
            %--------------------------------------------------------------
            A          = spm_clusters(xSwE.XYZ);
            mj          = find(A == A(i));
            xSwE.Z     = xSwE.Z(mj);
            xSwE.XYZ   = xSwE.XYZ(:,mj);
            xSwE.XYZmm = xSwE.XYZmm(:,mj);
        end

        %-Call 'list' functionality to produce table
        %------------------------------------------------------------------
        varargout = { swe_list('list',xSwE,hReg,Num,Dis,Str) };


    %======================================================================
    case 'txtlist'                                 %-Print ASCII text table
    %======================================================================
    % FORMAT swe_list('TxtList',TabDat,c)

        if nargin<2, error('Not enough input arguments.'); end
        if nargin<3, c = 1; else c = varargin{3}; end
        TabDat = varargin{2};

        %-Table Title
        %------------------------------------------------------------------
        fprintf('\n\nStatistics: %s\n',TabDat.tit)
        fprintf('%c',repmat('=',1,80)), fprintf('\n')

        %-Table header
        %------------------------------------------------------------------
        fprintf('%s\t',TabDat.hdr{1,c:end-1}), fprintf('%s\n',TabDat.hdr{1,end})
        fprintf('%s\t',TabDat.hdr{2,c:end-1}), fprintf('%s\n',TabDat.hdr{2,end})
        fprintf('%c',repmat('-',1,80)), fprintf('\n')

        %-Table data
        %------------------------------------------------------------------
        for i = 1:size(TabDat.dat,1)
            for mj=c:size(TabDat.dat,2)
                fprintf(TabDat.fmt{mj},TabDat.dat{i,mj});
                fprintf('\t')
            end
            fprintf('\n')
        end
        for i=1:max(1,11-size(TabDat.dat,1)), fprintf('\n'), end
        fprintf('%s\n',TabDat.str)
        fprintf('%c',repmat('-',1,80)), fprintf('\n')

        %-Table footer
        %------------------------------------------------------------------
        for i=1:size(TabDat.ftr,1)
            fprintf([TabDat.ftr{i,1} '\n'],TabDat.ftr{i,2});
        end
        fprintf('%c',repmat('=',1,80)), fprintf('\n\n')


    %======================================================================
    case 'xlslist'                                  %-Export table to Excel
    %======================================================================
    % FORMAT swe_list('XLSList',TabDat)

        if nargin<2, error('Not enough input arguments.'); end
        TabDat = varargin{2};

        d          = [TabDat.hdr(1:2,:);TabDat.dat];
        xyz        = d(3:end,end);
        xyz        = num2cell([xyz{:}]');
        d(:,end+1) = d(:,end);
        d(:,end+1) = d(:,end);
        d(3:end,end-2:end) = xyz;
        tmpfile    = [tempname '.xls'];
        xlswrite(tmpfile, d);
        winopen(tmpfile);

    %======================================================================
    case 'csvlist'            %-Export table to comma-separated values file
    %======================================================================
    % FORMAT swe_list('CSVList',TabDat)

        if nargin<2, error('Not enough input arguments.'); end
        TabDat = varargin{2};

        tmpfile = [tempname '.csv'];
        fid = fopen(tmpfile,'wt');
        fprintf(fid,[repmat('%s,',1,11) '%d,,\n'],TabDat.hdr{1,:});
        fprintf(fid,[repmat('%s,',1,11) '\n'],TabDat.hdr{2,:});
        fmt = TabDat.fmt;
        [fmt{2,:}] = deal(','); fmt = [fmt{:}];
        fmt(end:end+1) = '\n'; fmt = strrep(fmt,' ',',');
        for i=1:size(TabDat.dat,1)
            fprintf(fid,fmt,TabDat.dat{i,:});
        end
        fclose(fid);
        open(tmpfile);

    %======================================================================
    case 'setcoords'                                    %-Coordinate change
    %======================================================================
    % FORMAT swe_list('SetCoords',xyz,hAx,hReg)
        if nargin<3, error('Not enough input arguments.'); end
        hAx      = varargin{3};
        xyz      = varargin{2};
        UD       = get(hAx,'UserData');
        HlistXYZ = UD.HlistXYZ(ishandle(UD.HlistXYZ));

        %-Set all co-ord strings to black
        %------------------------------------------------------------------
        set(HlistXYZ,'Color','k');

        %-If co-ord matches a string, highlight it in red
        %------------------------------------------------------------------
        XYZ      = get(HlistXYZ,'UserData');
        if iscell(XYZ), XYZ = cat(2,XYZ{:}); end
        [tmp,i,d] = spm_XYZreg('NearestXYZ',xyz,XYZ);
        if d<eps
            set(HlistXYZ(i),'Color','r');
        end

    %======================================================================
    otherwise                                       %-Unknown action string
    %======================================================================
        error('Unknown action string')
end
