function [SwE,xSwE] = swe_getSPM(varargin)
% Compute specified and thresholded SwE parametric map for the SwE method.
% =========================================================================
% FORMAT [SwE,xSwE] = swe_getSPM;
% Query SwE in interactive mode.
%
% FORMAT [SwE,xSwE] = swe_getSPM(xSwE);
% -------------------------------------------------------------------------
%
% Query SwE in batch mode. See below for a description of fields that 
% may be present in xSwE input. Values for missing fields will be 
% queried interactively.
%
% xSwE       - structure containing spm, distribution & filtering 
%              details
% .swd       - SwE working directory - directory containing current 
%              SwE.mat
% .title     - title for comparison (string)
% .Z         - minimum of Statistics {filtered on u and k}
% .n         - conjunction number <= number of contrasts
% .STAT      - distribution {Z, T, X, F or P}
% .df        - degrees of freedom [df{interest}, df{residual}]
% .STATstr   - description string
% .Ic        - indices of contrasts (in SwE.xCon)
% .Im        - indices of masking contrasts (in xCon)
% .pm        - p-value for masking (uncorrected)
% .Ex        - flag for exclusive or inclusive masking
% .u         - height threshold
% .k         - extent threshold {voxels}
% .XYZ       - location of voxels {voxel coords}
% .XYZmm     - location of voxels {mm}
% .S         - search Volume {voxels}
% .R         - search Volume {resels}
% .FWHM      - smoothness {voxels}
% .M         - voxels -> mm matrix
% .iM        - mm -> voxels matrix
% .VOX       - voxel dimensions {mm} - column vector
% .DIM       - image dimensions {voxels} - column vector
% .Vspm      - Mapped statistic image(s)
% .Ps        - uncorrected P values in searched volume (for voxel FDR)
% .Pp        - uncorrected P values of peaks (for peak FDR)
% .Pc        - uncorrected P values of cluster extents (for cluster FDR)
% .uc        - 0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
% .thresDesc - description of height threshold (string)
%
% Required fields of SwE
%
% xVol   - structure containing details of volume analysed
%
% xX     - Design Matrix structure
%        - (see spm_spm.m for structure)
%
% xCon   - Contrast definitions structure array
%        - (see also spm_FcUtil.m for structure, rules & handling)
% .name  - Contrast name
% .STAT  - Statistic indicator character ('T', 'F' or 'P')
% .c     - Contrast weights (column vector contrasts)
% .X0    - Reduced design matrix data (spans design space under Ho)
%          Stored as coordinates in the orthogonal basis of xX.X from spm_sp
%          Extract using X0 = spm_FcUtil('X0',...
% .iX0   - Indicates how contrast was specified:
%          If by columns for reduced design matrix then iX0 contains the
%          column indices. Otherwise, it's a string containing the
%          spm_FcUtil 'Set' action: Usually one of {'c','c+','X0'}
% .X1o   - Remaining design space data (X1o is orthogonal to X0)
%          Stored as coordinates in the orthogonal basis of xX.X from spm_sp
%          Extract using X1o = spm_FcUtil('X1o',...
% .eidf  - Effective interest degrees of freedom (numerator df)
%        - Or effect-size threshold for Posterior probability
% .Vcon  - Name of contrast (for 'T's) or ESS (for 'F's) image
% .Vspm  - Name of SwE image
%
% Evaluated fields in xSwE (input)
%
% xSwE      - structure containing SwE, distribution & filtering details
% .swd      - SwE working directory - directory containing current SwE.mat
% .title    - title for comparison (string)
% .Ic       - indices of contrasts (in SwE.xCon)
% .n        - conjunction number <= number of contrasts
% .Im       - indices of masking contrasts (in xCon)
% .pm       - p-value for masking (uncorrected)
% .Ex       - flag for exclusive or inclusive masking
% .u        - height threshold
% .k        - extent threshold {voxels}
% .thresDesc - description of height threshold (string)
%
% In addition, the xCon structure is updated. For newly evaluated
% contrasts, SwE images (swe_vox_{T|F}stat_c{c#}) are written, along 
% with contrast (swe_vox_beta_c{c#}) images.
%
% For a parametric analysis the following is added to the xCon
% structure:
%
% .Vspm          - Name of SwE image
%
% For a non-parametric analysis the following are added to the xCon
% structure:
%
% .Vspm          - Name of SwE image
% .VspmFWEP      - Name of FWE P SwE image
% .VspmFDRP      - Name of FDR P SwE image
% .VspmUncP      - Name of Uncorrected P SwE image
% .VspmFWEP_clus - Name of FWE cluster P SwE image
% .Vedf          - Name of error degrees of freedom image
%
% The contrast images are the weighted sum of the parameter images,
% where the weights are the contrast weights, and are uniquely
% estimable since contrasts are checked for estimability by the
% contrast manager. These contrast images (for appropriate contrasts)
% are suitable summary images of an effect at this level, and can be
% used as input at a higher level when effecting a random effects
% analysis. (Note that the swe_vox_beta_c{c#} and
% swe_vox_{T|F}stat_c{c#} images are not suitable input for a higher
% level analysis.) See spm_RandFX.man for further details.
%
%__________________________________________________________________________
%
% swe_getSPM prompts for an SwE parametric map and applies thresholds {u & k}
% to a point list of voxel values (specified with their locations {XYZ})
% This allows the SwE map be displayed and characterized in terms of regionally
% significant effects by subsequent routines.
%
% For general linear model Y = XB + E with data Y, design matrix X,
% parameter vector B, and (independent) errors E, a contrast c'B of the
% parameters (with contrast weights c) is estimated by c'b, where b are
% the parameter estimates given by b=pinv(X)*Y.
%
% For a paramertic analysis, either single contrasts can be examined 
% or conjunctions of different contrasts. Contrasts are estimable 
% linear combinations of the parameters, and are specified using 
% the SwE contrast manager interface [swe_conman.m]. For a 
% non-parametric analysis, two contrasts are recorded; activation
% and deactivation for the contrast vector specified in the batch
% window. These are recorded a priori in a seperate function 
% with certain thresholds applied here [swe_contrasts_WB].
%
% SwE parametric maps are generated for the null hypotheses that 
% the contrast is zero (or zero vector in the case of F-contrasts). 
% See the help for the contrast manager [swe_conman.m] for a further 
% details on contrasts and contrast specification.
%
% A conjunction assesses the conjoint expression of multiple effects. The
% conjunction SwE is the minimum of the component SPMs defined by the
% multiple contrasts.  Inference on the minimum statistics can be
% performed in different ways.  Inference on the Conjunction Null (one or
% more of the effects null) is accomplished by assessing the minimum as
% if it were a single statistic; one rejects the conjunction null in
% favor of the alternative that k=nc, that the number of active effects k
% is equal to the number of contrasts nc.  No assumptions are needed on
% the dependence between the tests.
%
% Another approach is to make inference on the Global Null (all effects
% null).  Rejecting the Global Null of no (u=0) effects real implies an
% alternative that k>0, that one or more effects are real.   A third
% Intermediate approach, is to use a null hypothesis of no more than u
% effects are real.  Rejecting the intermediate null that k<=u implies an
% alternative that k>u, that more than u of the effects are real.
%
% The Global and Intermediate nulls use results for minimum fields which
% require the SPMs to be identically distributed and independent. Thus,
% all component SwE maps must be either SwE{t}'s, or SwE{F}'s with the same
% degrees of freedom. Independence is roughly guaranteed for large
% degrees of freedom (and independent data) by ensuring that the
% contrasts are "orthogonal". Note that it is *not* the contrast weight
% vectors per se that are required to be orthogonal, but the subspaces of
% the data space implied by the null hypotheses defined by the contrasts
% (c'pinv(X)). Furthermore, this assumes that the errors are
% i.i.d. (i.e. the estimates are maximum likelihood or Gauss-Markov. This
% is the default in spm_spm).
%
% To ensure approximate independence of the component SwE maps in the case of
% the global or intermediate null, non-orthogonal contrasts are serially
% orthogonalised in the order specified, possibly generating new
% contrasts, such that the second is orthogonal to the first, the third
% to the first two, and so on.  Note that significant inference on the
% global null only allows one to conclude that one or more of the effects
% are real.  Significant inference on the conjunction null allows one to
% conclude that all of the effects are real.
%
% Masking simply eliminates voxels from the current contrast if they
% do not survive an uncorrected p value (based on height) in one or
% more further contrasts.  No account is taken of this masking in the
% statistical inference pertaining to the masked contrast.
%
% The SwE map is subject to thresholding on the basis of height (u) and the
% number of voxels comprising its clusters {k}. The height threshold is
% specified as above in terms of an [un]corrected p value or
% statistic.  Clusters can also be thresholded on the basis of their
% spatial extent. If you want to see all voxels simply enter 0.  In this
% instance the 'set-level' inference can be considered an 'omnibus test'
% based on the number of clusters that obtain.
%
% BAYESIAN INFERENCE AND PPMS - POSTERIOR PROBABILITY MAPS
%
% If conditional estimates are available (and your contrast is a T
% contrast) then you are asked whether the inference should be 'Bayesian'
% or 'classical' (using GRF).  If you choose Bayesian the contrasts are of
% conditional (i.e. MAP) estimators and the inference image is a
% posterior probability map (PPM).  PPMs encode the probability that the
% contrast exceeds a specified threshold.  This threshold is stored in
% the xCon.eidf.  Subsequent plotting and tables will use the conditional
% estimates and associated posterior or conditional probabilities.
%
% see swe_results_ui.m for further details of the SwE results section.
% see also swe_contrasts.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Modified version of spm_getSPM
% Written by B. Guillaume
% Version Info:  $Format:%ci$ $Format:%h$

%-GUI setup
%--------------------------------------------------------------------------
%spm_help('!ContextHelp',mfilename)
spm('Pointer','Arrow')

%-Select SwE.mat & note SwE results directory
%--------------------------------------------------------------------------
if nargin
    xSwE = varargin{1};
end
try
    swd = xSwE.swd;
catch
    swd = '.';
end

%-Preliminaries...
%==========================================================================

%-Load SwE.mat
%--------------------------------------------------------------------------
try
    load(fullfile(swd,'SwE.mat'));
catch
    error(['Cannot read ' fullfile(swd,'SwE.mat')]);
end

SwE.swd = swd;

%-Change directory so that relative filenames are valid
%--------------------------------------------------------------------------
cd(SwE.swd);

%-Check the model has been estimated
%--------------------------------------------------------------------------
try
    SwE.xVol.XYZ;
catch
    
    %-Check the model has been estimated
    %----------------------------------------------------------------------
    error( 'This model has not been estimated.');
    
end

% check format of data
file_ext = swe_get_file_extension(SwE.xY.P{1});
isMat    = strcmpi(file_ext,'.mat');
isCifti  = strcmpi(file_ext,'.dtseries.nii') ||  strcmpi(file_ext,'.dtscalar.nii');

if ~isMat && ~isCifti
  isMeshData = spm_mesh_detect(SwE.xY.VY);
end

xX   = SwE.xX;                      %-Design definition structure
XYZ  = SwE.xVol.XYZ;                %-XYZ coordinates
S    = SwE.xVol.S;                  %-search Volume {voxels}
% R    = SwE.xVol.R;                  %-search Volume {resels}
if isMat
  M = SwE.xVol.M;
  VOX = [];
  clear xSwE
else
  M    = SwE.xVol.M(1:3,1:3);         %-voxels to mm matrix
  VOX  = sqrt(diag(M'*M))';           %-voxel dimensions
end

% Tolerance for comparing real numbers for WB analyses
% Use a value < to the smallest WB p-value as it will be used to include WB p-values equal to alpha
if isfield(SwE, 'WB') 
  tol = 0.1 / (SwE.WB.nB + 1);
end

% check the data and other files have valid filenames
%--------------------------------------------------------------------------
%something here occurs and the paths to spm and swe toolboxes disappear????
% try, SwE.xY.VY       = spm_check_filename(SwE.xY.VY);       end
% try, SwE.xVol.VRpv   = spm_check_filename(SwE.xVol.VRpv);   end
% try, SwE.Vbeta       = spm_check_filename(SwE.Vbeta);       end
% try, SwE.Vcov_vis    = spm_check_filename(SwE.Vcov_vis);    end
% try, SwE.Vcov_beta   = spm_check_filename(SwE.Vcov_beta);   end
% try, SwE.Vcov_beta_g = spm_check_filename(SwE.Vcov_beta_g); end %here seems the problem
% try, SwE.VM          = spm_check_filename(SwE.VM);          end

%==========================================================================
% - C O N T R A S T S ,   S P M    C O M P U T A T I O N ,    M A S K I N G
%==========================================================================

%-Get contrasts
%--------------------------------------------------------------------------
try
  xCon = SwE.xCon;
  % check if the Uncorrected p-value image is correctly set to the non-parametric version for WB (for retro-compatibility)
  if isfield(SwE, 'WB') && ~exist('OCTAVE_VERSION','builtin') && ~contains(xCon(1).VspmUncP.fname, '-WB')
    for i = 1:numel(xCon)
      SwE.xCon(i).VspmUncP = spm_vol(sprintf('swe_vox_%cstat_lp%s_c%.2d%s', SwE.WB.stat, '-WB', i, file_ext));
    end
    % save the modified SwE.mat
    if spm_check_version('matlab','7') >=0
      save('SwE.mat', 'SwE', '-V6');
    else
      save('SwE.mat', 'SwE');
    end
    xCon = SwE.xCon; 
  end
catch
  if isfield(SwE, 'WB') && ~exist('OCTAVE_VERSION','builtin')
    SwE = swe_contrasts_WB(SwE);
    % save SwE with xCon appended to it. This is important for future call of swe_getSPM for a specific Ic
    if spm_check_version('matlab','7') >=0
      save('SwE.mat', 'SwE', '-V6');
    else
      save('SwE.mat', 'SwE');
    end
    xCon = SwE.xCon;    
  else
    xCon = {};
  end
end

try
    Ic        = xSwE.Ic;
catch
    
    % If we're not doing wild bootstrap and not in octave, ask for a contrast.
    if ~isfield(SwE, 'WB') && ~exist('OCTAVE_VERSION','builtin')
        [Ic,xCon] = swe_conman(SwE,'T&F',Inf,...
                               '    Select contrasts...',' for conjunction',1);
    % If we're in octave, assume we already have a contrast.
    elseif exist('OCTAVE_VERSION','builtin')
        Ic = 1;
    % If we're doing WB, we already have a contrast. We just need to record it.
    else
        if numel(xCon) == 2
            Ic = spm_input('Contrast Type','+1','b','Activation|Deactivation',[1,2],1);
        else
            Ic = 1;
        end
    end
end
if isempty(xCon)
    % figure out whether new contrasts were defined, but not selected
    % do this by comparing length of SwE.xCon to xCon, remember added
    % indices to run spm_contrasts on them as well
    try
        noxCon = numel(SwE.xCon);
    catch
        noxCon = 0;
    end
    IcAdd = (noxCon+1):numel(xCon);
else
    IcAdd = [];
end

nc        = length(Ic);  % Number of contrasts

%-Allow user to extend the null hypothesis for conjunctions
%
% n: conjunction number
% u: Null hyp is k<=u effects real; Alt hyp is k>u effects real
%    (NB Here u is from Friston et al 2004 paper, not statistic thresh).
%                  u         n
% Conjunction Null nc-1      1     |    u = nc-n
% Intermediate     1..nc-2   nc-u  |    #effects under null <= u
% Global Null      0         nc    |    #effects under alt  > u,  >= u+1
%----------------------------------+---------------------------------------
if nc > 1
    try
        n = xSwE.n;
    catch
        if nc==2
            But = 'Conjunction|Global';      Val=[1 nc];
        else
            But = 'Conj''n|Intermed|Global'; Val=[1 NaN nc];
        end
        n = spm_input('Null hyp. to assess?','+1','b',But,Val,1);
        if isnan(n)
            if nc == 3
                n = nc - 1;
            else
                n = nc - spm_input('Effects under null ','0','n1','1',nc-1);
            end
        end
    end
else
    n = 1; 
end

% not sure we want to do that with the SwE (commented for now)
%-Enforce orthogonality of multiple contrasts for conjunction
% (Orthogonality within subspace spanned by contrasts)
%--------------------------------------------------------------------------
% if nc > 1 && n > 1 && ~spm_FcUtil('|_?',xCon(Ic), xX.xKXs)
%     
%     OrthWarn = 0;
%     
%     %-Successively orthogonalise
%     %-NB: This loop is peculiarly controlled to account for the
%     %     possibility that Ic may shrink if some contrasts disappear
%     %     on orthogonalisation (i.e. if there are colinearities)
%     %----------------------------------------------------------------------
%     i = 1;
%     while(i < nc), i = i + 1;
%         
%         %-Orthogonalise (subspace spanned by) contrast i w.r.t. previous
%         %------------------------------------------------------------------
%         oxCon = spm_FcUtil('|_',xCon(Ic(i)), xX.xKXs, xCon(Ic(1:i-1)));
%         
%         %-See if this orthogonalised contrast has already been entered
%         % or is colinear with a previous one. Define a new contrast if
%         % neither is the case.
%         %------------------------------------------------------------------
%         d     = spm_FcUtil('In',oxCon,xX.xKXs,xCon);
%         
%         if spm_FcUtil('0|[]',oxCon,xX.xKXs)
%             
%             %-Contrast was colinear with a previous one - drop it
%             %--------------------------------------------------------------
%             Ic(i) = [];
%             i     = i - 1;
%             
%         elseif any(d)
%             
%             %-Contrast unchanged or already defined - note index
%             %--------------------------------------------------------------
%             Ic(i) = min(d);
%             
%         else
%             
%             %-Define orthogonalised contrast as new contrast
%             %--------------------------------------------------------------
%             OrthWarn   = OrthWarn + 1;
%             conlst     = sprintf('%d,',Ic(1:i-1));
%             oxCon.name = sprintf('%s (orth. w.r.t {%s})', xCon(Ic(i)).name,...
%                                   conlst(1:end-1));
%             xCon       = [xCon, oxCon];
%             Ic(i)      = length(xCon);
%         end
%         
%     end % while...
%     
%     if OrthWarn
%         warning('SwE:ConChange','%d contrasts orthogonalized',OrthWarn)
%     end
%     
%     SwE.xCon = xCon;
% end % if nc>1...
SwE.xCon = xCon;

%-Apply masking
%--------------------------------------------------------------------------
try
  Mask = ~isempty(xSwE.Im) * (isnumeric(xSwE.Im) + 2*iscellstr(xSwE.Im));
catch
  % Mask = spm_input('mask with other contrast(s)','+1','y/n',[1,0],2);
  if isMat
    % no masking optionfor mat format
    Mask = 0;
  elseif isfield(SwE, 'WB')
    % for now, the post-hoc masking is disabled for the WB
    % It may be added later.
    Mask = 0;
    % Mask = spm_input('apply masking','+1','b','none|image',[0,2],1);
  else
    Mask = spm_input('apply masking','+1','b','none|contrast|image',[0,1,2],1);
  end
end
if Mask == 1

  %-Get contrasts for masking
  %----------------------------------------------------------------------
  try
    Im = xSwE.Im;
  catch
    [Im,xCon] = swe_conman(SwE,'T&F',-Inf,...
      'Select contrasts for masking...',' for masking',1);
  end

  %-Threshold for mask (uncorrected p-value)
  %----------------------------------------------------------------------
  try
   pm = xSwE.pm;
  catch
    pm = spm_input('uncorrected mask p-value','+1','r',0.05,1,[0,1]);
  end

  %-Inclusive or exclusive masking
  %----------------------------------------------------------------------
  try
    Ex = xSwE.Ex;
  catch
    Ex = spm_input('nature of mask','+1','b','inclusive|exclusive',[0,1],1);
  end

elseif Mask == 2

  %-Get mask images
  %----------------------------------------------------------------------
  try
    Im = xSwE.Im;
  catch
    if isMat
      Im = cellstr(spm_select([1 Inf],'mat','Select mask image(s)'));
    else
      [Im, sts] = spm_select([1 Inf],{'image','mesh'},'Select mask image(s)');
      if ~sts, Im = []; else Im = cellstr(Im); end
    end
  end

  %-Inclusive or exclusive masking
  %----------------------------------------------------------------------
  try
    Ex = xSwE.Ex;
  catch
    Ex = spm_input('nature of mask','+1','b','inclusive|exclusive',[0,1],1);
  end

  pm = [];

else
  Im = [];
  pm = [];
  Ex = [];
end

%-Create/Get title string for comparison
%--------------------------------------------------------------------------
if isMat
  titlestr = xCon(Ic).name;
else
  if nc == 1
    str  = xCon(Ic).name;
  else
    str  = [sprintf('contrasts {%d',Ic(1)),sprintf(',%d',Ic(2:end)),'}'];
    if n == nc
      str = [str ' (global null)'];
    elseif n == 1
      str = [str ' (conj. null)'];
    else
      str = [str sprintf(' (Ha: k>=%d)',(nc-n)+1)];
    end
  end
  if Ex
    mstr = 'masked [excl.] by';
  else
    mstr = 'masked [incl.] by';
  end
  if isnumeric(Im)
    if length(Im) == 1
      str = sprintf('%s (%s %s at p=%g)',str,mstr,xCon(Im).name,pm);
    elseif ~isempty(Im)
      str = [sprintf('%s (%s {%d',str,mstr,Im(1)),...
        sprintf(',%d',Im(2:end)),...
        sprintf('} at p=%g)',pm)];
    end
  elseif iscellstr(Im) && numel(Im) > 0
    [pf,nf,ef] = spm_fileparts(Im{1});
    str  = sprintf('%s (%s %s',str,mstr,[nf ef]);
    for i=2:numel(Im)
      [pf,nf,ef] = spm_fileparts(Im{i});
      str =[str sprintf(', %s',[nf ef])];
    end
    str = [str ')'];
  end
end

if ~isMat
  try
    titlestr = xSwE.title;
    if isempty(titlestr)
        titlestr = str;
    end
  catch
    titlestr = spm_input('title for comparison','+1','s',str);
  end
end

if ~isMat
    % Ask whether to do additional voxelwise or clusterwise inference.
    try
        infType = xSwE.infType;
    catch 
        if isfield(SwE, 'WB')
            if isfield(SwE.WB, 'TFCE')
                infType = spm_input('inference type','+1','b','voxelwise|clusterwise|TFCE',[0,1,2],3);
            else
                infType = spm_input('inference type','+1','b','voxelwise|clusterwise',[0,1],2);
            end
        else
            infType = spm_input('inference type','+1','b','voxelwise|clusterwise',[0,1],1);
        end
    end
    
    if isfield(SwE, 'WB')
        % Work out the original form of inference performed. This will tell us
        % which maps have already been generated. Most importantly, whether we
        % can do FWE p value clusterwise inference.
        if SwE.WB.voxelWise
            orig_infType = 'vox';
        elseif SwE.WB.clusterWise 
            orig_infType = 'clus';
        else
            orig_infType = 'tfce';
        end
    end
    
end

%-Compute & store contrast parameters, contrast/ESS images, & SwE images
%==========================================================================
SwE.xCon = xCon;
alreadyComputed = all(~cellfun(@isempty,{xCon(Ic).Vspm}));

if isnumeric(Im)
    SwE  = swe_contrasts(SwE, unique([Ic, Im, IcAdd]));
else
    SwE  = swe_contrasts(SwE, unique([Ic, IcAdd]));
end
xCon     = SwE.xCon;
STAT     = xCon(Ic(1)).STAT;
VspmSv   = cat(1,xCon(Ic).Vspm);

    
    
%-Check conjunctions - Must be same STAT w/ same df
%--------------------------------------------------------------------------
if (nc > 1) && (any(diff(double(cat(1,xCon(Ic).STAT)))) || ...
                any(abs(diff(cat(1,xCon(Ic).eidf))) > 1))
    error('illegal conjunction: can only conjoin SPMs of same STAT & df');
end


%-Degrees of Freedom and STAT string describing marginal distribution
%--------------------------------------------------------------------------
%dFWHM=SwE.xVol.FWHM * SwE.xVol.M(1:3,1:3);
%df=[xCon(Ic(1)).eidf, (SwE.Subj.nSubj-SwE.xX.pB)*((1+2*(SwE.vFWHM(1)/dFWHM(1))^2)*(1+2*(SwE.vFWHM(2)/dFWHM(2))^2)*(1+2*(SwE.vFWHM(3)/dFWHM(3))^2))^0.5-xCon(Ic(1)).eidf+1]; %df=[xCon(Ic(1)).eidf nSubj-pB-xCon(Ic(1)).eidf-1]; %df     = [xCon(Ic(1)).eidf xX.erdf];
if nc > 1
    if n > 1
        str = sprintf('^{%d \\{Ha:k\\geq%d\\}}',nc,(nc-n)+1);
    else
        str = sprintf('^{%d \\{Ha:k=%d\\}}',nc,(nc-n)+1);
    end
else
    str = '';
end

% We display the equivalent statistics.
switch STAT
    case 'T' 
        STATstr = sprintf('%c','Z',str);
    case 'F'
        STATstr = sprintf('%c','X',str);
end
%-Compute (unfiltered) spm pointlist for masked conjunction requested
%==========================================================================
fprintf('\t%-32s: %30s','SPM computation','...initialising')            %-#


%-Compute conjunction as minimum of SPMs
%--------------------------------------------------------------------------
Z     = Inf;
for i = Ic
  if isMat
    load(xCon(i).Vspm);
    Z = min(Z,equivalentScore);
    clear equivalentScore
  else
    Z = min(Z, swe_data_read(xCon(i).Vspm, 'xyz', XYZ));
  end
      
end


%-Copy of Z and XYZ before masking, for later use with FDR
%--------------------------------------------------------------------------
XYZum = XYZ;
Zum   = Z;

%-Compute mask and eliminate masked voxels
%--------------------------------------------------------------------------
for i = 1:numel(Im)

    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...masking')           %-#
    if isnumeric(Im)
      if isMat
        load(xCon(Im(i)).Vspm);
        Mask = equivalentScore;
        clear equivalentScore
      else
        Mask = swe_data_read(xCon(Im(i)).Vspm, 'xyz', XYZ);
      end
      switch xCon(Im(i)).STAT
        case 'T'
          um   = swe_invNcdf(1-pm);
        case 'F'
          um   = spm_invXcdf(1-pm,1);
      end
      if Ex
        Q = Mask <= um;
      else
        Q = Mask >  um;
      end
    else
      if isMat
        Mask = importdata(Im{i});
      else
        v = swe_data_hdr_read(Im{i});
        Mask = swe_data_read(v, 'xyz', v.mat\SwE.xVol.M*[XYZ; ones(1,size(XYZ,2))]);
      end
      Q = Mask ~= 0 & ~isnan(Mask);
      if Ex, Q = ~Q; end
    end
    if ~isMat
      XYZ   = XYZ(:,Q);
    end
    Z     = Z(Q);
    if isempty(Q)
        fprintf('\n')                                                   %-#
        warning('SwE:NoVoxels','No voxels survive masking at p=%4.2f',pm);
        break
    end
end


  %==========================================================================
  % - H E I G H T   &   E X T E N T   T H R E S H O L D S
  %==========================================================================
if ~isMat
  u   = -Inf;        % height threshold
  k   = 0;           % extent threshold {voxels}
  clustWise = 'None';% Type of clusterwise inference to be performed

  if  spm_mesh_detect(xCon(Ic(1)).Vspm)
    G = export(gifti(SwE.xVol.G),'patch');
  end

  %-Height threshold - classical inference
  %--------------------------------------------------------------------------
  if STAT ~= 'P'
      
      % Get the equivalent statistic
      switch STAT
          
          case 'T'
              
              eSTAT = 'Z';
              
          case 'F'
              
              eSTAT = 'X';
              
      end
      
      % If we are doing voxelwise inference on a parametric.
      if ~isfield(SwE, 'WB') && infType == 0
          
          %-Get height threshold
          %----------------------------------------------------------------
          fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...height threshold')  %-#
          try
              thresDesc = xSwE.thresDesc;
          catch
              % For non WB we only have FDR.
              str = 'FDR|none';

              thresDesc = spm_input('p value adjustment to control','+1','b',str,[],1);
          end
          
          switch thresDesc

              case 'FDR' % False discovery rate
                  %--------------------------------------------------------
                  try
                      u = xSwE.u;
                  catch
                      u = spm_input('p value (FDR)','+0','r',0.05,1,[0,1]);
                  end
                  thresDesc = ['p<' num2str(u) ' (' thresDesc ')'];
                  switch STAT
                      case 'T'
                         u = spm_uc_FDR(u,Inf,'Z',n,VspmSv,0); 
                      case 'F'
                         u = spm_uc_FDR(u,[1 1],'X',n,VspmSv,0); 
                  end

              case 'none'  % No adjustment: p for conjunctions is p of the conjunction SwE
                  %--------------------------------------------------------
                  try
                      u = xSwE.u;
                  catch
                      u = spm_input(['threshold {',eSTAT,' or p value}'],'+0','r',0.001,1);
                  end
                  if u <= 1
                      thresDesc = ['p<' num2str(u) ' (unc.)'];
                      switch STAT
                          case 'T'
                              u  = swe_invNcdf(1-u^(1/n));
                          case 'F'
                              u  = spm_invXcdf(1-u^(1/n),1);
                      end
                  else
                      '';
                  end

              otherwise
                  %--------------------------------------------------------------
                  fprintf('\n');                                              %-#
                  error('Unknown control method "%s".',thresDesc);
                  

          end % switch thresDesc
          
          %-Compute p-values for topological and voxel-wise FDR (all search voxels)
          %----------------------------------------------------------------
          fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...for voxelFDR')  %-#
          switch STAT
              case 'T'
                  Ps = (1-spm_Ncdf(Zum)).^n; 
              case 'F'
                  Ps = (1-spm_Xcdf(Zum,1)).^n;
          end
          
          up  = NaN;
          Pp  = NaN;
          uc  = NaN;
          ue  = NaN;
          Pc  = [];
          uu = [];
          
          Q      = find(Z > u);
      
      % If we are doing clusterwise inference on a parametric.
      elseif ~isfield(SwE, 'WB') && infType == 1
          
          % Record what type of clusterwise inference we are doing.
          clustWise = 'Uncorr';
          
          % No adjustment: p for conjunctions is p of the conjunction SwE
          %--------------------------------------------------------------
          try
              u = xSwE.u;
          catch
              u = spm_input(['threshold {',eSTAT,' or p value}'],'+1','r',0.001,1);
          end
          if u <= 1
              thresDesc = ['p<' num2str(u) ' (unc.)'];
              switch STAT
                  case 'T'
                      u  = swe_invNcdf(1-u^(1/n));
                  case 'F'
                      u  = spm_invXcdf(1-u^(1/n),1);
              end
          else
              thresDesc = '';
          end
          
          %-Compute p-values for topological and voxel-wise FDR (all search voxels)
          %----------------------------------------------------------------
          fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...for voxelFDR')  %-#
          switch STAT
              case 'T'
                  Ps = (1-spm_Ncdf(Zum)).^n; 
              case 'F'
                  Ps = (1-spm_Xcdf(Zum,1)).^n;
          end
          
          up  = NaN;
          Pp  = NaN;
          uc  = NaN;
          ue  = NaN;
          Pc  = [];
          uu = [];
          
          Q      = find(Z > u);

      % If we are doing voxelwise inference on a WB.
      elseif isfield(SwE, 'WB') && infType == 0
          
          %-Get height threshold
          %----------------------------------------------------------------------
          fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...height threshold')  %-#
          try
              thresDesc = xSwE.thresDesc;
          catch
              % For WB we have FWE or FDR.
              str = 'FWE|FDR|none';
              thresDesc = spm_input('p value adjustment to control','+1','b',str,[],1);
          end
          
          switch thresDesc

              case 'FWE' % Family-wise false positive rate
                  % This is performed on the voxelwise FWE P value map
                  %--------------------------------------------------------
                  try
                      pu = xSwE.u;
                  catch
                      pu = spm_input('p value (FWE)','+0','r',0.05,1,[0,1]);
                  end
                  thresDesc = ['p<=' num2str(pu) ' (' thresDesc ')'];
                  
                  FWE_ps = 10.^-swe_data_read(xCon(Ic).VspmFWEP,'xyz', XYZ);
                  
                  % When thresholding on WB FWER p-values, we should include those = to pu
                  % Here, we are using a - tol < b instead of a <= b due to numerical errors
                  % tol was set to 0.1/(nB+1) in order to make sure it is smaller than the smallest WB p-value
                  Q = find(FWE_ps - tol  < pu);
                  
                  % Obtain the exclusive statistic threshold. This will be the (1-pu)th
                  % percentile of the max. statistic distribution
                  if Ic == 1
                    maxScore = sort(SwE.WB.maxScore);
                  elseif Ic == 2
                    maxScore = sort(-SwE.WB.minScore);
                  else
                    error("Unknown contrast");
                  end
                  u = maxScore( ceil( (1-pu) * (SwE.WB.nB+1) ) );

              case 'FDR' % False discovery rate
                  % This is performed on the FDR P value map
                  %--------------------------------------------------------
                  try
                      pu = xSwE.u;
                  catch
                      pu = spm_input('p value (FDR)','+0','r',0.05,1,[0,1]);
                  end
                  thresDesc = ['p<=' num2str(pu) ' (' thresDesc ')'];
                  
                  % select the WB FDR p-values within the mask
                  FDR_ps = 10.^-swe_data_read(xCon(Ic).VspmFDRP, 'xyz', XYZ);

                  % Here, a parametric score threshold u would differ from voxel to voxel
                  % Thus, setting it to NaN
                  u = NaN
                  
                  % inclusive thresholding for WB
                  Q = find(FDR_ps - tol < pu);

              case 'none'  % No adjustment: p for conjunctions is p of the conjunction SwE
                  % This should be performed on the uncorrected WB p-values
                  %--------------------------------------------------------
                  try
                      pu = xSwE.u;
                  catch
                      pu = spm_input(['threshold {p value}'],'+0','r',0.001,1,[0,1]);
                  end
                  thresDesc = ['p<=' num2str(pu) ' (unc.)'];
                  % select the WB unc. p-values within the mask
                  unc_ps = 10.^-swe_data_read(xCon(Ic).VspmUncP, 'xyz', XYZ);

                  % Here, a parametric score threshold u would differ from voxel to voxel
                  % Thus, setting it to NaN
                  u = NaN
                  
                  % inclusive thresholding for WB
                  Q = find(unc_ps - tol < pu);

              otherwise
                  %--------------------------------------------------------------
                  fprintf('\n');                                              %-#
                  error('Unknown control method "%s".',thresDesc);
          end
          
          up  = NaN;
          Pp  = NaN;
          uc  = NaN;
          ue  = NaN;
          Pc  = [];
          uu = [];

      % If we are doing clusterwise WB.
      elseif isfield(SwE, 'WB') && infType == 1
          
          fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...height threshold')  %-#
          try
              thresDesc = xSwE.thresDesc;
          catch
              if strcmp(orig_infType, 'clus')
                % For WB we have FWE or p/k specification if the user
                % originally ran a clusterwise analysis.
                str = 'FWE|none';
                thresDesc = spm_input('extent p value adjustment to control','+1','b',str,[],1);
              else
                % If the user originally ran voxelwise or TFCE we have no
                % FWE clusterwise map from the bootstrap.
                thresDesc = 'none';
              end
          end
          
          switch thresDesc
              
              case 'FWE' % Family-wise false positive rate
                  % This is performed on the FWE P value map
                  %--------------------------------------------------------
                  
                  % Record what type of clusterwise inference we are doing.
                  clustWise = 'FWE';
                  
                  % For WB we have performed our thresholding and worked out of 
                  % regions of activation. 
                  if Ic == 1
                    locActVox = SwE.WB.clusterInfo.LocActivatedVoxels;
                  elseif Ic == 2
                    locActVox = SwE.WB.clusterInfo.LocActivatedVoxelsNeg;
                  else
                    error("Unknown contrast");
                  end
                  [~, index]=ismember(XYZ',locActVox','rows');
                  Q=find(index~=0);
                  
                  % We also should record the cluster forming threshold that was
                  % used.
                  pu = SwE.WB.clusterInfo.primaryThreshold;
                  thresDesc = ['p<' num2str(pu) ' (unc.)'];
                  if strcmp(STAT, 'T')
                      u = norminv(1-pu);
                  else
                      u = chi2inv(1-pu, 1);
                  end
                  
                  % We should display the cluster forming threshold that
                  % was used.
                  spm_input(['threshold {p value}'],...
                      '+1','b',['(pre-set: P=' num2str(pu) ')'],[0],0)
                  
              case 'none'  % No adjustment: p for conjunctions is p of the conjunction SwE
                  % This should be performed on the uncorrected WB p-values
                  %--------------------------------------------------------
                  % Record what type of clusterwise inference we are doing.
                  clustWise = 'Uncorr';
                  
                  % Cluster-forming threshold.
                  try
                      pu = xSwE.u;
                  catch
                      pu = spm_input(['threshold {p value}'],'+0','r',0.001,1,[0,1]);
                  end
                  thresDesc = ['p<=' num2str(pu) ' (unc.)'];
                  % select the WB unc. p-values within the mask
                  unc_ps = 10.^-swe_data_read(xCon(Ic).VspmUncP, 'xyz', XYZ);
  
                  % Here, a parametric score threshold u would differ from voxel to voxel
                  % Thus, setting it to NaN
                  u = NaN
                  
                  % inclusive thresholding for WB
                  Q = find(unc_ps - tol < pu);
                  
                  up  = NaN;
                  Pp  = NaN;
                  uc  = NaN;
                  ue  = NaN;
                  Pc  = [];
                  uu = [];
                            
          end
          
      % If we are doing TFCE.    
      else
          
          % Remember we are not doing clusterwise.
          clustWise = 'None';
          
          % Ask user for TFCE FWE alpha.
          pt = spm_input('p value (TFCE FWE)','+0','r',0.05,1,[0,1]);
          thresDesc = ['p<=' num2str(pt) ' (FWE)'];
          
          % Get Tfce Fwe P-values.
          % In older version of the toolbox, the max TFCE scores were not saved.
          % Thus, to avoid retro-compatibility issues, we cannot threshold using 
          % the (1-pt)th percentile of the max distribution, but only using the FWER p-values
          tfp = 10.^-swe_data_read(xCon(Ic).VspmTFCEFWEP, 'xyz', XYZ);
          
          up  = NaN;
          Pp  = NaN;
          uc  = NaN;
          ue  = NaN;
          Pc  = [];
          uu = [];

          % When thresholding on WB p-values, we should include those = to pu
          % Here, we are using a - tol < b instead of a <= b due to numerical errors
          % tol was set to 0.1/(nB+1) in order to make sure it is smaller than the smallest WB p-value
          Q = find(tfp - tol < pt);
          
      end

  end

  %-Apply height threshold
  %--------------------------------------------------------------------------
  Z      = Z(:,Q);
  if ~isMat
    XYZ    = XYZ(:,Q);
  end
	if isMat
		strDataType = 'data elements';
	else
		if spm_mesh_detect(xCon(Ic(1)).Vspm)
			strDataType = 'vertices'
		else 
			strDataType = 'voxels'; 
		end
	end
  if isempty(Q)
      fprintf('\n');
      warning('SwE:NoVoxels','No %s survive thresholding', strDataType);
  end
  
  % If we are doing clusterwise ask for threshold.
  if ~strcmp(clustWise, 'None')
      %-Extent threshold (disallowed for conjunctions)
      %--------------------------------------------------------------------------
      if ~isempty(XYZ) && nc == 1 

          fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...extent threshold')  %-#
          
          % Uncorrected threshold.
          if strcmp(clustWise, 'Uncorr')
              %-Get extent threshold [default = 0]
              %----------------------------------------------------------------------
              try
                  k = xSwE.k;
              catch
                  k = spm_input('& extent threshold {voxels}','+1','r',0,1,[0,Inf]);
              end

              %-Calculate extent threshold filtering
							%----------------------------------------------------------------------
              if isCifti
                A = swe_cifti_clusters(SwE.cifti, XYZ(1,:));              
              elseif  ~spm_mesh_detect(xCon(Ic(1)).Vspm)
                A = spm_clusters(XYZ);
              else
                T = false(SwE.xVol.DIM');
                T(XYZ(1,:)) = true;
                A = spm_mesh_clusters(G,T)';
                A = A(XYZ(1,:));
              end
              Q     = [];
              for i = 1:max(A)
                  j = find(A == i);
                  if length(j) >= k; Q = [Q j]; end
              end
          end
          
          % FWE WB threshold.
          if strcmp(clustWise, 'FWE')
              %-Get extent threshold [default = 0]
              %----------------------------------------------------------------------
              try
                  fwep_c = xSwE.fwep_c;
              catch
                  fwep_c = spm_input('& extent threshold {FWE P}','+1','r',0.05,1,[0,1]);
              end

              %-Calculate extent threshold filtering
              %----------------------------------------------------------------------
              % recompute the clusters as they may have been reduced due to post-hoc masking
              if isCifti
                A = swe_cifti_clusters(SwE.cifti, XYZ(1,:));              
              elseif  ~spm_mesh_detect(xCon(Ic(1)).Vspm)
                A = spm_clusters(XYZ);
              else
                T = false(SwE.xVol.DIM');
                T(XYZ(1,:)) = true;
                A = spm_mesh_clusters(G,T)';
                A = A(XYZ(1,:));
            	end
              clusIndices = unique(A);
              
              % recompute the p-values as they might have increased due to post-hoc masking
              if Ic == 1
                for i = 1:length(clusIndices)
                  ind = ( A==clusIndices(i) );
                  ps_fwe(ind) = sum( SwE.WB.clusterInfo.maxClusterSize >= sum(ind) ) / (SwE.WB.nB + 1);
                end
              elseif Ic == 2
                for i = 1:length(clusIndices)
                  ind = ( A==clusIndices(i) );
                  ps_fwe(ind) = sum( SwE.WB.clusterInfo.maxClusterSizeNeg >= sum(ind) ) / (SwE.WB.nB + 1);
                end
              else
                  error("unknown contrast");
              end

              % select only the voxels surviving the FWER threshold 
              % Here, we use an inclusive p-value threshold to be consistent with an exclusive cluster threshold
              Q = find(ps_fwe - tol < fwep_c);
              
              % The exclusive threshold k should be the (1-fwep_c)th percentile of the max cluster size distribution
              if Ic == 1
                maxClusterSize = sort(SwE.WB.clusterInfo.maxClusterSize);
              elseif Ic == 2
                maxClusterSize = sort(SwE.WB.clusterInfo.maxClusterSizeNeg);
              else
                error("Unknown contrast");
              end
              k = maxClusterSize( ceil( (1-fwep_c) * (SwE.WB.nB+1) ) );

          end

          % ...eliminate voxels
          %----------------------------------------------------------------------
          Z     = Z(:,Q);
          XYZ   = XYZ(:,Q);
          if isempty(Q)
              fprintf('\n');                                                  %-#
              warning('SwE:NoVoxels','No %s survive cluster extent threshoding', strDataType);
          end
          
      else

          k = 0;

      end 
  end
end 

% Doftype
if isfield(SwE.type,'modified')
    dof_type = SwE.type.modified.dof_mo;
else
    dof_type = SwE.type.classic.dof_cl;        
end

%==========================================================================
% - E N D
%==========================================================================
if alreadyComputed && isMat
  fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...already done previously')                %-#
else
  fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#
end
spm('Pointer','Arrow')

%-Assemble output structures of unfiltered data
%==========================================================================
if isMat
  XYZmm     = [];
  Z         = [];
  u         = [];
  k         = [];
  thresDesc = [];
else
  XYZmm = SwE.xVol.M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];
end
xSwE   = struct( ...
            'swd',        swd,...
            'title',      titlestr,...
            'Z',          Z,...
            'n',          n,...
            'STAT',       STAT,...
            'STATstr',    STATstr,...
            'Ic',         Ic,...
            'Im',         {Im},...
            'pm',         pm,...
            'Ex',         Ex,...
            'u',          u,...
            'k',          k,...
            'XYZ',        XYZ,...
            'XYZmm',      XYZmm,...
            'S',          SwE.xVol.S,...
            'M',          SwE.xVol.M,...
            'iM',         SwE.xVol.iM,...
            'DIM',        SwE.xVol.DIM,...
            'VOX',        VOX,...
            'Vspm',       VspmSv,...
            'thresDesc',  thresDesc,...
            'WB',         0,...
            'dofType',    dof_type,...
            'nPredict',   size(SwE.xX.X, 2),...
            'df_Con',     rank(xCon(Ic).c),...
            'SS',         SwE.SS);
        
if isfield(SwE.type, 'modified')
    xSwE.nSubj_g    = SwE.Gr.nSubj_g;
    xSwE.max_nVis_g = SwE.Vis.max_nVis_g;
    xSwE.min_nVis_g = SwE.Vis.min_nVis_g;
end

if dof_type == 0
  xSwE.edf = xCon(Ic).edf;
else
  xSwE.Vedf = cat(1,xCon(Ic).Vedf);
end

% For WB analyses we have already computed uncorrected, FDR, FWE and
% cluster-FWE P values at this point.
if isfield(SwE, 'WB')
    
    % Bootstrap variables.
    xSwE.WB = 1;
    xSwE.nB = SwE.WB.nB;
    
    % Volumes
    xSwE.VspmUncP = cat(1,xCon(Ic).VspmUncP);
    xSwE.VspmFDRP = cat(1,xCon(Ic).VspmFDRP);
    xSwE.VspmFWEP = cat(1,xCon(Ic).VspmFWEP);
    if SwE.WB.clusterWise
        xSwE.VspmFWEP_clus = cat(1,xCon(Ic).VspmFWEP_clus);
    end
    if isfield(SwE.WB, 'TFCE')
        xSwE.VspmTFCE = cat(1,xCon(Ic).VspmTFCE);
        xSwE.VspmTFCEP = cat(1,xCon(Ic).VspmTFCEP);
        xSwE.VspmTFCEFWEP = cat(1,xCon(Ic).VspmTFCEFWEP);
        xSwE.TFCEanaly = 1;
        xSwE.TFCE = SwE.WB.TFCE;
    else
        xSwE.TFCEanaly = 0;
    end
    if infType == 2
        xSwE.TFCEthresh = 1;
    else
        xSwE.TFCEthresh = 0;
    end
    
    % Uncorrected P values.
    Ps_vol = swe_data_hdr_read(xSwE.VspmUncP);
    Ps = swe_data_read(Ps_vol);
    Ps = 10.^(-Ps(~isnan(Ps)));
    xSwE.Ps = Ps;
    
    % 95% percentiles
    if Ic == 1
      maxScore = sort(SwE.WB.maxScore);
    elseif Ic == 2
      maxScore = sort(-SwE.WB.minScore);
    else
      error("Unknown contrast");
    end
    xSwE.Pfv = maxScore(ceil(0.95*(xSwE.nB+1))); % Voxelwise FWE P 
    if SwE.WB.clusterWise
      if Ic == 1
        maxClusterSize = sort(SwE.WB.clusterInfo.maxClusterSize);
      elseif Ic == 2
        maxClusterSize = sort(SwE.WB.clusterInfo.maxClusterSizeNeg);
      else
        error("Unknown contrast");
      end
      xSwE.Pfc = maxClusterSize(ceil(0.95*(xSwE.nB+1))); % Clusterwise FWE P      
    end
    
    % edf
    xSwE.Vedf       = cat(1,xCon(Ic).Vedf);
    
end

% Record clusterwise FWE P value if there is one.
if ~isMat
    xSwE.infType = infType;
    xSwE.clustWise = clustWise;
    if strcmp(clustWise, 'FWE')
        xSwE.fwep_c = fwep_c;
    end
end
 %             'R',        SwE.xVol.R,...
%             'FWHM',     SwE.xVol.FWHM,...

xSwE.isCifti = isCifti;
if isCifti
  xSwE.cifti = SwE.cifti;
  tmp = 0;
  if numel(SwE.cifti.surfaces) > 0
    for i = 1:numel(SwE.cifti.surfaces)
      indInSurface = SwE.cifti.surfaces{i}.off + (1:numel(SwE.cifti.surfaces{i}.iV));
      isSurviving = ismember(indInSurface, SwE.xVol.XYZ(1,:));
      tmp = tmp + sum(isSurviving);
    end
  end
  xSwE.S_surf = tmp;
  if numel(SwE.cifti.volume) > 0
    isSurviving = ismember(SwE.cifti.volume.indices, SwE.xVol.XYZ(1,:));
    xSwE.S_vol = sum(isSurviving);
  else
    xSwE.S_vol = 0;
  end
end
          
%-RESELS per voxel (density) if it exists
%--------------------------------------------------------------------------
try, xSwE.VRpv = SwE.xVol.VRpv; end
try
    xSwE.units = SwE.xVol.units;
catch
    try, xSwE.units = varargin{1}.units; end
end

%-Topology for surface-based inference
%--------------------------------------------------------------------------
if spm_mesh_detect(xCon(Ic(1)).Vspm)
  xSwE.G     = G;
  xSwE.XYZmm = xSwE.G.vertices(xSwE.XYZ(1,:),:)';
end

%-p-values for topological and voxel-wise FDR
%--------------------------------------------------------------------------
try, xSwE.Ps    = Ps;             end  % voxel   FDR
try, xSwE.Pp    = Pp;             end  % peak    FDR
try, xSwE.Pc    = Pc;             end  % cluster FDR

%-0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
%--------------------------------------------------------------------------
try, xSwE.uc    = [uu up ue uc];  end
