function [SwE,xSwE] = swe_getSPM(varargin)
% Compute a specified and thresholded SwE for the SwE method
% FORMAT [SwE,xSwE] = swe_getSPM;
% Query SwE in interactive mode.
%
% FORMAT [SwE,xSwE] = swe_getSPM(xSwE);
% Query SwE in batch mode. See below for a description of fields that may
% be present in xSwE input. Values for missing fields will be queried
% interactively.
%
% xSwE      - structure containing spm, distribution & filtering details
% .swd      - SwE working directory - directory containing current SwE.mat
% .title    - title for comparison (string)
% .Z        - minimum of Statistics {filtered on u and k}
% .n        - conjunction number <= number of contrasts
% .STAT     - distribution {Z, T, X, F or P}
% .df       - degrees of freedom [df{interest}, df{residual}]
% .STATstr  - description string
% .Ic       - indices of contrasts (in SwE.xCon)
% .Im       - indices of masking contrasts (in xCon)
% .pm       - p-value for masking (uncorrected)
% .Ex       - flag for exclusive or inclusive masking
% .u        - height threshold
% .k        - extent threshold {voxels}
% .XYZ      - location of voxels {voxel coords}
% .XYZmm    - location of voxels {mm}
% .S        - search Volume {voxels}
% .R        - search Volume {resels}
% .FWHM     - smoothness {voxels}
% .M        - voxels -> mm matrix
% .iM       - mm -> voxels matrix
% .VOX      - voxel dimensions {mm} - column vector
% .DIM      - image dimensions {voxels} - column vector
% .Vspm     - Mapped statistic image(s)
% .Ps       - uncorrected P values in searched volume (for voxel FDR)
% .Pp       - uncorrected P values of peaks (for peak FDR)
% .Pc       - uncorrected P values of cluster extents (for cluster FDR)
% .uc       - 0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
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
% contrasts, SwE images (spmT_????) are written, along with
% contrast (con_????) images for SwE{T}'s, or Extra
% Sum-of-Squares images (ess_????}) for SwE{F}'s.
%
% The contrast images are the weighted sum of the parameter images,
% where the weights are the contrast weights, and are uniquely
% estimable since contrasts are checked for estimability by the
% contrast manager. These contrast images (for appropriate contrasts)
% are suitable summary images of an effect at this level, and can be
% used as input at a higher level when effecting a random effects
% analysis. (Note that the ess_????.{img,hdr} and
% SwE{T,F}_????.{img,hdr} images are not suitable input for a higher
% level analysis.) See spm_RandFX.man for further details.
%
%__________________________________________________________________________
%
% swe_getSPM prompts for an SwE and applies thresholds {u & k}
% to a point list of voxel values (specified with their locations {XYZ})
% This allows the SwE be displayed and characterized in terms of regionally
% significant effects by subsequent routines.
%
% For general linear model Y = XB + E with data Y, design matrix X,
% parameter vector B, and (independent) errors E, a contrast c'B of the
% parameters (with contrast weights c) is estimated by c'b, where b are
% the parameter estimates given by b=pinv(X)*Y.
%
% Either single contrasts can be examined or conjunctions of different
% contrasts. Contrasts are estimable linear combinations of the
% parameters, and are specified using the SwE contrast manager
% interface [swe_conman.m]. SPMs are generated for the null hypotheses
% that the contrast is zero (or zero vector in the case of
% F-contrasts). See the help for the contrast manager [swe_conman.m]
% for a further details on contrasts and contrast specification.
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
% all component SPMs must be either SwE{t}'s, or SwE{F}'s with the same
% degrees of freedom. Independence is roughly guaranteed for large
% degrees of freedom (and independent data) by ensuring that the
% contrasts are "orthogonal". Note that it is *not* the contrast weight
% vectors per se that are required to be orthogonal, but the subspaces of
% the data space implied by the null hypotheses defined by the contrasts
% (c'pinv(X)). Furthermore, this assumes that the errors are
% i.i.d. (i.e. the estimates are maximum likelihood or Gauss-Markov. This
% is the default in spm_spm).
%
% To ensure approximate independence of the component SPMs in the case of
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
% The SwE is subject to thresholding on the basis of height (u) and the
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
% see spm_results_ui.m for further details of the SwE results section.
% see also spm_contrasts.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%Modified version of spm_getSPM
%Written by B. Guillaume

%-GUI setup
%--------------------------------------------------------------------------
%spm_help('!ContextHelp',mfilename)
if ~exist('OCTAVE_VERSION','builtin')
  spm('Pointer','Arrow')
end

%-Select SwE.mat & note SwE results directory
%--------------------------------------------------------------------------
if nargin
    xSwE = varargin{1};
end
try
    swd = xSwE.swd;
    sts = 1;
catch
    if ~exist('OCTAVE_VERSION','builtin')
        [spmmatfile, sts] = spm_select(1,'^SwE\.mat$','Select SwE.mat');
        swd = spm_str_manip(spmmatfile,'H');
    else
        sts = 1;
        swd = '.';
    end
end
if ~sts, SwE = []; xSwE = []; return; end

%-Preliminaries...
%==========================================================================

%-Load SwE.mat
%--------------------------------------------------------------------------
try
    load(fullfile(swd,'SwE.mat'));
catch
    error(['Cannot read ' fullfile(swd,'SwE.mat')]);
end

if isfield(SwE, 'WB')
  msg = {'No result display feature is currently availabe for the Wild Bootstrap. The results are written into -log10(p-values) images (e.g., see lP_FWE+.img for FWE-corrected p-values) into the folder used for the analysis. You can display and threshold these images using SPM or an alternative software package.'
  '\n\nBelow, a quick description of all possible outputs:'
  '\n - swe_vox_lp-WB_c0001{neg}: Images of -log10(voxel-wise uncorrected non-parametric P-values, positive or negative.)'
  '\n - swe_vox_lpFWE-WB_c0001{neg}: Images of -log10(voxel-wise FWE-corrected non-parametric P-values, positive or negative). Here, FWE-corrected non-parametric P-values are the proportion of the wild bootstrap distribution for the maximal statistic which exceeds the statistic image at the voxel.'
  '\n - swe_vox_lpFDR-WB_c0001{neg}: Images of -log10(voxel-wise FDR-corrected non-parametric P-values, positive or negative). They are computed based by applying BH-fdr correction on lP+ & lP-.'
  '\n - swe_clus_lpFWE-WB_c0001{neg}: Images of -log10(cluster-wise FWE-corrected non-parametric P-values, positive or negative). Note that, the -log10(p-values) of each formed cluster is repeated at each voxel belonging to this cluster. More iinformation about the formed cluters can be found in the "SwE.mat" file in the field SwE.WB.clusterInfo.'
  '%s'};
  error(strcat(msg{1:6}),'');
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
[~,~,file_ext] = fileparts(SwE.xY.P{1});
isMat          = strcmpi(file_ext,'.mat');

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
try xCon = SwE.xCon; catch, xCon = {}; end

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
        xCon = SwE.xCon;
    % If we're doing WB, we already have a contrast. We just need to record it.
    else
        Ic = 1;
        xCon = struct('name', ['swe_', SwE.WB.stat, 'stat-WB_01'],...
                      'STAT', SwE.WB.stat,...
                      'c', SwE.WB.con);
        if SwE.WB.stat == 'T' 
            if SwE.WB.clusterWise == 1
                xCon.Vspm = struct('fname', '');
            end 
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
      Mask = 0;
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
        Im = cellstr(spm_select([1 Inf],'image','Select mask image(s)'));
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
  try
      titlestr = xSwE.title;
      if isempty(titlestr)
          titlestr = str;
      end
  catch
      titlestr = spm_input('title for comparison','+1','s',str);
  end
end

%-Compute & store contrast parameters, contrast/ESS images, & SwE images
%==========================================================================
SwE.xCon = xCon;
alreadyComputed = ~isempty(xCon(Ic).Vspm);

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
df     = [1 1];
if nc > 1
    if n > 1
        str = sprintf('^{%d \\{Ha:k\\geq%d\\}}',nc,(nc-n)+1);
    else
        str = sprintf('^{%d \\{Ha:k=%d\\}}',nc,(nc-n)+1);
    end
else
    str = '';
end

switch STAT
    case 'T' 
        STATstr = sprintf('%c','T',str);
    case 'F'
        STATstr = sprintf('%c','F',str);
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
    Z = min(Z,spm_get_data(xCon(i).Vspm,XYZ));
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
        Mask = spm_get_data(xCon(Im(i)).Vspm,XYZ);
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
        v = spm_vol(Im{i});
        Mask = spm_get_data(v,v.mat\SwE.xVol.M*[XYZ; ones(1,size(XYZ,2))]);
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


  %-Height threshold - classical inference
  %--------------------------------------------------------------------------
  if STAT ~= 'P'

      %-Get height threshold
      %----------------------------------------------------------------------
      fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...height threshold')  %-#
      try
          thresDesc = xSwE.thresDesc;
      catch
            str = 'FDR|none';
  %         if topoFDR
  %             str = 'FWE|none';
  %         else
  %             str = 'FWE|FDR|none';
  %         end
          thresDesc = spm_input('p value adjustment to control','+1','b',str,[],1);
      end

      switch thresDesc

  %         case 'FWE' % Family-wise false positive rate
  %             %--------------------------------------------------------------
  %             try
  %                 u = xSwE.u;
  %             catch
  %                 u = spm_input('p value (FWE)','+0','r',0.05,1,[0,1]);
  %             end
  %             thresDesc = ['p<' num2str(u) ' (' thresDesc ')'];
  %             switch STAT
  %                 case 'T'
  %                     u   = spm_uc(u,[0,0],'Z',R,n,S);
  %                 case 'F'
  %                     u   = spm_uc(u,[0,1],'X',R,n,S);
  %             end


          case 'FDR' % False discovery rate
              %--------------------------------------------------------------
  %             if topoFDR,
  %                 fprintf('\n');                                          %-#
  %                 error('Change defaults.stats.topoFDR to use voxel FDR');
  %             end
              try
                  u = xSwE.u;
              catch
                  u = spm_input('p value (FDR)','+0','r',0.05,1,[0,1]);
              end
              thresDesc = ['p<' num2str(u) ' (' thresDesc ')'];
              switch STAT
                  case 'T'
                     u = spm_uc_FDR(u,df,'Z',n,VspmSv,0); 
                  case 'F'
                     u = spm_uc_FDR(u,df,'X',n,VspmSv,0); 
              end

          case 'none'  % No adjustment: p for conjunctions is p of the conjunction SwE
              %--------------------------------------------------------------
              try
                  u = xSwE.u;
              catch
                  u = spm_input(['threshold {',STAT,' or p value}'],'+0','r',0.001,1);
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
                  thresDesc = [STAT '=' num2str(u) ];
              end



          otherwise
              %--------------------------------------------------------------
              fprintf('\n');                                              %-#
              error('Unknown control method "%s".',thresDesc);

      end % switch thresDesc
      %-Compute p-values for topological and voxel-wise FDR (all search voxels)
      %----------------------------------------------------------------------
      fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...for voxelFDR')  %-#
      switch STAT
          case 'T'
              Ps = (1-spm_Ncdf(Zum)).^n; 
          case 'F'
              Ps = (1-spm_Xcdf(Zum,df(2))).^n;
      end
      %-Peak FDR
      %----------------------------------------------------------------------
  %     switch STAT
  %         case 'T'
  %             [up,Pp] = spm_uc_peakFDR(0.05,df,'Z',R,n,Zum,XYZum,u);
  %         case 'F'
  %             [up,Pp] = spm_uc_peakFDR(0.05,df,'X',R,n,Zum,XYZum,u);
  %     end
          up  = NaN;
          Pp  = NaN;
      %-Cluster FDR
      %----------------------------------------------------------------------
  %     if STAT == 'T' && n == 1
  %         V2R        = 1/prod(SwE.xVol.FWHM(SwE.xVol.DIM > 1));
  %         [uc,Pc,ue] = spm_uc_clusterFDR(0.05,df,STAT,R,n,Zum,XYZum,V2R,u);
  %     else
          uc  = NaN;
          ue  = NaN;
          Pc  = [];
  %     end

      %-Peak FWE
      %----------------------------------------------------------------------
  %     switch STAT
  %         case 'T'
  %             uu      = spm_uc(0.05,[0 0],'Z',R,n,S);
  %         case 'F'
  %             uu      = spm_uc(0.05,[0 1],'X',R,n,S);
  %     end
      uu = [];

  end

  %-Calculate height threshold filtering
  %--------------------------------------------------------------------------
  Q      = find(Z > u);

  %-Apply height threshold
  %--------------------------------------------------------------------------
  Z      = Z(:,Q);
  if ~isMat
    XYZ    = XYZ(:,Q);
  end
  if isempty(Q)
      fprintf('\n');                                                      %-#
      warning('SwE:NoVoxels','No voxels survive masking at p=%4.2f',pm);
  end


  %-Extent threshold (disallowed for conjunctions)
  %--------------------------------------------------------------------------
  if ~isempty(XYZ) && nc == 1

      fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...extent threshold')  %-#

      %-Get extent threshold [default = 0]
      %----------------------------------------------------------------------
      try
          k = xSwE.k;
      catch
          k = spm_input('& extent threshold {voxels}','+1','r',0,1,[0,Inf]);
      end

      %-Calculate extent threshold filtering
      %----------------------------------------------------------------------
      A     = spm_clusters(XYZ);
      Q     = [];
      for i = 1:max(A)
          j = find(A == i);
          if length(j) >= k; Q = [Q j]; end
      end

      % ...eliminate voxels
      %----------------------------------------------------------------------
      Z     = Z(:,Q);
      XYZ   = XYZ(:,Q);
      if isempty(Q)
          fprintf('\n');                                                  %-#
          warning('SwE:NoVoxels','No voxels survive masking at p=%4.2f',pm);
      end

  else

      k = 0;

  end % (if ~isempty(XYZ))
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
            'swd',      swd,...
            'title',    titlestr,...
            'Z',        Z,...
            'n',        n,...
            'STAT',     STAT,...
            'df',       df,...
            'STATstr',  STATstr,...
            'Ic',       Ic,...
            'Im',       {Im},...
            'pm',       pm,...
            'Ex',       Ex,...
            'u',        u,...
            'k',        k,...
            'XYZ',      XYZ,...
            'XYZmm',    XYZmm,...
            'S',        SwE.xVol.S,...
            'M',        SwE.xVol.M,...
            'iM',       SwE.xVol.iM,...
            'DIM',      SwE.xVol.DIM,...
            'VOX',      VOX,...
            'Vspm',     VspmSv,...
            'thresDesc',thresDesc);

 %             'R',        SwE.xVol.R,...
%             'FWHM',     SwE.xVol.FWHM,...
          
          
%-RESELS per voxel (density) if it exists
%--------------------------------------------------------------------------
try, xSwE.VRpv = SwE.xVol.VRpv; end
try
    xSwE.units = SwE.xVol.units;
catch
    try, xSwE.units = varargin{1}.units; end
end

%-p-values for topological and voxel-wise FDR
%--------------------------------------------------------------------------
try, xSwE.Ps    = Ps;             end  % voxel   FDR
try, xSwE.Pp    = Pp;             end  % peak    FDR
try, xSwE.Pc    = Pc;             end  % cluster FDR

%-0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
%--------------------------------------------------------------------------
try, xSwE.uc    = [uu up ue uc];  end
