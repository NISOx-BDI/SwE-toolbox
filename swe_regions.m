function [Y,xY] = swe_regions(xSwE,SwE,hReg,xY)
% VOI time-series extraction of adjusted data (& local eigenimage analysis)
% FORMAT [Y,xY] = swe_regions(xSwE,SwE,hReg,[xY])
%
% xSwE   - structure containing specific SwE, distribution & filtering details
% SwE    - structure containing generic analysis details
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Y      - first scaled eigenvariate of VOI {i.e. weighted mean}
% xY     - VOI structure
%       xY.xyz          - centre of VOI {mm}
%       xY.name         - name of VOI
%       xY.Ic           - contrast used to adjust data (0 - no adjustment)
%       xY.Sess         - session index
%       xY.def          - VOI definition
%       xY.spec         - VOI definition parameters
%       xY.str          - VOI description as a string
%       xY.XYZmm        - Co-ordinates of VOI voxels {mm}
%       xY.y            - [whitened and filtered] voxel-wise data
%       xY.u            - first eigenvariate {scaled - c.f. mean response}
%       xY.v            - first eigenimage
%       xY.s            - eigenvalues
%       xY.X0           - [whitened] confounds (including drift terms)
%
% Y and xY are also saved in VOI_*.mat in the SwE working directory.
% (See spm_getSwE for details on the SwE & xSwE structures)
%
% FORMAT [Y,xY] = swe_regions('Display',[xY])
%
% xY     - VOI structure or filename
%
%__________________________________________________________________________
%
% swe_regions extracts a representative time course from voxel data in
% terms of the first eigenvariate of the filtered and adjusted response in
% all suprathreshold voxels within a specified VOI centred on the current
% MIP cursor location. Responses are adjusted by removing variance that
% can be predicted by the null space of the F contrast specified (usually
% an F-contrast testing for all effects of interest).
%
% If temporal filtering has been specified, then the data will be filtered.
% Similarly for whitening. Adjustment is with respect to the null space of
% a selected contrast, or can be omitted.
%
% For a VOI of radius 0, the [adjusted] voxel time-series is returned, and
% scaled to have a 2-norm of 1. The actual [adjusted] voxel time series can
% be extracted from xY.y, and will be the same as the [adjusted] data
% returned by the plotting routine (spm_graph.m) for the same contrast.
%__________________________________________________________________________
% Copyright (C) 1999-2016 Wellcome Trust Centre for Neuroimaging
% Karl Friston

% Code adapted from swe_regions by Tom Maullin (28/11/2018)
% Version Info:  $Format:%ci$ $Format:%h$
% Based on: swe_regions.m 6923 2016-11-04 15:35:12Z guillaume

% Add the design matrix into the xKXs field to match that of an SwE.mat
SwE.xX.xKXs.X = SwE.xX.X;

%-Shortcut for VOI display
%--------------------------------------------------------------------------
if nargin && ischar(xSwE) && strcmpi(xSwE,'display')
    if nargin > 1
        xY = SwE;
    else
        [xY, sts] = spm_select(1,'^VOI.*\.mat$');
        if ~sts, return; end
    end
    if ischar(xY), load(xY); else Y = xY.u; end
    display_VOI(xY);
    return;
end

if nargin < 4, xY = []; end

%-Get figure handles
%--------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), noGraph = 1; else noGraph = 0; end
header = get(Finter,'Name');
set(Finter,'Name','VOI time-series extraction');
if ~noGraph, spm_figure('GetWin','Graphics'); end

%-Find nearest voxel [Euclidean distance] in point list
%--------------------------------------------------------------------------
if isempty(xSwE.XYZmm)
    spm('alert!','No suprathreshold voxels!',mfilename,0);
    Y = []; xY = [];
    return
end
try
    xyz    = xY.xyz;
catch
    xyz    = spm_XYZreg('NearestXYZ',...
             spm_XYZreg('GetCoords',hReg),xSwE.XYZmm);
    xY.xyz = xyz;
end

% and update GUI location
%--------------------------------------------------------------------------
if spm_mesh_detect(SwE.xY.VY)
    spm_XYZreg('SetCoords',xyz,hReg,1);
else
    spm_XYZreg('SetCoords',xyz,hReg);
end


%-Get adjustment options and VOI name
%--------------------------------------------------------------------------
if ~noGraph
    if ~isempty(xY.xyz)
        posstr = sprintf('at [%3.0f %3.0f %3.0f]',xY.xyz);
    else
        posstr = '';
    end
    spm_input(posstr,1,'d','VOI time-series extraction');
end

if ~isfield(xY,'name')
    xY.name    = spm_input('name of region','!+1','s','VOI');
end

if ~isfield(xY,'Ic')
    q(1)   = 0;
    Con    = {'<don''t adjust>'};
    q(2)   = NaN;
    Con{2} = '<adjust for everything>';
    for i = 1:length(SwE.xCon)
        if strcmp(SwE.xCon(i).STAT,'F')
            q(end + 1) = i;
            Con{end + 1} = SwE.xCon(i).name;
        end
    end
    if numel(Con) == 2
        warning('No F-contrast has been defined: are you sure?');
    end
    i     = spm_input('adjust data for (select contrast)','!+1','m',Con);
    xY.Ic = q(i);
end

%-If fMRI data then ask user to select session
%--------------------------------------------------------------------------
if isfield(SwE,'Sess') && ~isfield(xY,'Sess')
    s       = length(SwE.Sess);
    if s > 1
        s   = spm_input('which session','!+1','n1',s,Inf);
    end
    xY.Sess = s;
end

%-Specify VOI
%--------------------------------------------------------------------------
xY.M = xSwE.M;
[xY, xY.XYZmm, Q] = spm_ROI(xY, xSwE.XYZmm);
try, xY = rmfield(xY,'M'); end
try, xY = rmfield(xY,'rej'); end

if isempty(xY.XYZmm)
    warning('Empty region.');
    [Y, xY.y, xY.u, xY.v, xY.s, xY.X0] = deal([]);
    return;
end

%-Perform time-series extraction to all sessions if Inf is entered
%--------------------------------------------------------------------------
if isfield(SwE,'Sess') && isfield(xY,'Sess') && isinf(xY.Sess)
    if length(SwE.Sess) == 1
        xY.Sess = 1;
    else
        for i=1:length(SwE.Sess)
            xY.Sess = i;
            [tY{i},txY(i)] = swe_regions(xSwE,SwE,hReg,xY);
        end
        Y = tY; xY = txY;
        return;
    end
end

%-Extract required data from results files
%==========================================================================
spm('Pointer','Watch')

%-Get raw data, whiten and filter
%--------------------------------------------------------------------------
y        = swe_data_read(SwE.xY.VY,'xyz',xSwE.XYZ(:,Q));
%y        = spm_filter(SwE.xX.K,SwE.xX.W*y);


%-Computation
%==========================================================================

%-Remove null space of contrast
%--------------------------------------------------------------------------
if xY.Ic ~= 0

    %-Parameter estimates: beta = xX.pKX*xX.K*y
    %----------------------------------------------------------------------
    beta  = swe_data_read(SwE.Vbeta,'xyz',xSwE.XYZ(:,Q));

    %-subtract Y0 = XO*beta,  Y = Yc + Y0 + e
    %----------------------------------------------------------------------
    if ~isnan(xY.Ic)
        y = y - spm_FcUtil('Y0',SwE.xCon(xY.Ic),SwE.xX.xKXs.X,beta);
    else
        y = y - SwE.xX.xKXs.X * beta;
    end

end

%-Confounds
%--------------------------------------------------------------------------
xY.X0     = SwE.xX.xKXs.X(:,[SwE.xX.iB SwE.xX.iG]);

%-Extract session-specific rows from data and confounds
%--------------------------------------------------------------------------
try
    i     = SwE.Sess(xY.Sess).row;
    y     = y(i,:);
    xY.X0 = xY.X0(i,:);
end

% and add session-specific filter confounds
%--------------------------------------------------------------------------
try
    xY.X0 = [xY.X0 SwE.xX.K(xY.Sess).X0];
end
try
    xY.X0 = [xY.X0 SwE.xX.K(xY.Sess).KH]; % Compatibility check
end

%-Remove null space of X0
%--------------------------------------------------------------------------
xY.X0     = xY.X0(:,any(xY.X0));


%-Compute regional response in terms of first eigenvariate
%--------------------------------------------------------------------------
[m,n]   = size(y);
if m > n
    [v,s,v] = svd(y'*y);
    s       = diag(s);
    v       = v(:,1);
    u       = y*v/sqrt(s(1));
else
    [u,s,u] = svd(y*y');
    s       = diag(s);
    u       = u(:,1);
    v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u*sqrt(s(1)/n);

%-Set in structure
%--------------------------------------------------------------------------
xY.y    = y;
xY.u    = Y;
xY.v    = v;
xY.s    = s;

%-Display VOI weighting and eigenvariate
%==========================================================================
if ~noGraph
    if isfield(SwE.xY,'RT')
        display_VOI(xY, SwE.xY.RT);
    else
        display_VOI(xY);
    end
end

%-Save
%==========================================================================
str = ['VOI_' xY.name '.mat'];
if isfield(xY,'Sess') && isfield(SwE,'Sess')
    str = sprintf('VOI_%s_%i.mat',xY.name,xY.Sess);
end
save(fullfile(SwE.swd,str),'Y','xY', spm_get_defaults('mat.format'))

cmd = 'swe_regions(''display'',''%s'')';
fprintf('   VOI saved as %s\n',spm_file(fullfile(SwE.swd,str),'link',cmd));

%-Reset title
%--------------------------------------------------------------------------
set(Finter,'Name',header);
spm('Pointer','Arrow')


%==========================================================================
% function display_VOI(xY,TR)
%==========================================================================
function display_VOI(xY,TR)
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);
figure(Fgraph);
fullsize = isempty(get(Fgraph,'children'));

%-Show position
%--------------------------------------------------------------------------
if fullsize, subplot(2,1,1); else subplot(2,2,3); end
spm_dcm_display(xY)
if fullsize, title(['Region: ' xY.name]); end

%-Show dynamics
%--------------------------------------------------------------------------
if fullsize, subplot(2,1,2); else subplot(2,2,4); end
if nargin == 2
    plot(TR*[1:length(xY.u)],xY.u);
    str = 'time \{seconds\}';
else
    plot(xY.u);
    str = 'scan';
end
title(['1st eigenvariate: ' xY.name],'FontSize',10)
if strcmpi(xY.def,'mask')
    posstr = sprintf('from mask %s', spm_file(xY.spec.fname,'filename'));
else
    posstr = sprintf('at [%3.0f %3.0f %3.0f]',xY.xyz);
end
str = { str;' ';...
    sprintf('%d voxels in VOI %s',size(xY.y,2),posstr);...
    sprintf('Variance: %0.2f%%',100*xY.s(1)/sum(xY.s))};
xlabel(str)
axis tight square
