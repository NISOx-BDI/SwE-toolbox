function [Y,y,beta,Bcov] = swe_graph(xSwE,SwE,hReg)
% Graphical display of adjusted data
% =========================================================================
% FORMAT [Y y beta Bcov] = swe_graph(xSwE,SPM,hReg)
% -------------------------------------------------------------------------
% Inputs: 
%  - xSwE   - structure containing SPM, distributional & filtering details
%             about the excursion set
%  - SwE    - structure containing generic details about the analysis
%  - hReg   - handle of MIP register or [x y z] coordinates
% -------------------------------------------------------------------------
% Outputs:
%  - Y      - fitted   data for the selected voxel
%  - y      - adjusted data for the selected voxel
%  - beta   - parameter estimates (ML or MAP)
%  - Bcov   - Covariance of parameter estimates (ML or conditional)
%
% See swe_getSPM for details.
%__________________________________________________________________________
%
% swe_graph is a Callback script that uses the structures above to:  (i)
% send adjusted (y) and fitted data (Y), for the selected voxel, to the
% workspace and (ii) provide graphics for:
%
% a) Contrasts of parameter estimates (e.g. activations) and their
% standard error.
%
% b) Fitted and adjusted responses that can be plotted against time, scan,
% or an indicator variable in the design matrix.
%
% c) (fMRI only).  Evoked responses using the basis functions to give
% impulse responses that would have been seen in the absence of other
% effects. The PSTH (peristimulus-time histogram) option provides a finite
% impulse response (FIR) estimate of the trial-specific evoked response as
% a function of peristimulus time.  This is estimated by refitting a
% convolution model to the selected voxel using an FIR basis set.  This is
% simply a set of small boxes covering successive time bins after trial
% onset.  The width of each bin is usually the TR.  This option provides a
% more time-resolved quantitative characterisation of the evoked
% hemodynamic response.  However, it should not be over-interpreted because
% inference is usually made using a simpler and more efficient basis set
% (e.g., canonical hrf, or canonical plus time derivative).
%
% Getting adjusted data:
% Ensuring the data are adjusted properly can be important (e.g. in
% constructing explanatory variables such as in a psychophysiological
% interaction). To remove or correct for specific effects, specify an
% appropriate F contrast and simply plot the fitted (and adjusted)
% responses after selecting that F contrast. The vectors Y (fitted) and y
% (adjusted) in the workspace will now be corrected for the effects in the
% reduced design matrix (X0) specified in the contrast manager with the
% column indices (iX0) of the confounds in this adjustment.
%
% Plotting data:
% All data and graphics use filtered/whitened data and residuals. In PET
% studies the parameter estimates and the fitted data are often the same
% because the explanatory variables are simply indicator variables taking
% the value of one.  Only contrasts previously defined can be plotted. This
% ensures that the parameters plotted are meaningful even when there is
% collinearity among the design matrix subpartitions.
%
% Selecting contrasts used for PPMs will automatically give plots
% based on conditonal estimates.
%
% The structure     contrast.contrast      = cbeta;
%                   contrast.standarderror = SE;
%                   contrast.interval      = 2*CI;
%
% is assigned in base workspace for plots of contrasts and their error.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Modified version of swe_graph
% Modified by Bryan Guillaume
% Version Info:  $Format:%ci$ $Format:%h$


%-Get Graphics figure handle
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');


%-Delete previous axis and their pagination controls (if any)
%--------------------------------------------------------------------------
swe_results_ui('Clear',Fgraph,2);


%-Find nearest voxel [Euclidean distance] in point list & update GUI
%--------------------------------------------------------------------------
if isempty(xSwE.XYZmm)
    spm('alert!','No suprathreshold voxels!',mfilename,0);

    Y = []; y = []; beta = []; Bcov = [];
    return
end

if numel(hReg) == 1
    xyz = spm_XYZreg('GetCoords',hReg);
else
    xyz = hReg;
end
[xyz,i] = spm_XYZreg('NearestXYZ',xyz,xSwE.XYZmm);
if numel(hReg) == 1, spm_XYZreg('SetCoords',xyz,hReg); end
XYZ     = xSwE.XYZ(:,i);


%-Find out what to plot
%==========================================================================
Cplot = {   'Contrast estimates and 95% C.I.',...
            'Fitted responses',...
            'Event-related responses',...
            'Parametric responses',...
            'Volterra Kernels'};


% ensure options are appropriate
%--------------------------------------------------------------------------
try
    Sess  = SwE.Sess;
catch
    Cplot = Cplot(1);
end
Cplot  = Cplot{spm_input('Plot',-1,'m',Cplot)};

switch Cplot

    % select contrast if
    %----------------------------------------------------------------------
    case {'Contrast estimates and 95% C.I.','Fitted responses'}

        % determine which contrast
        %------------------------------------------------------------------
        ind_tcon = [];
        for i = 1:size(SwE.xCon,2)
            if size(SwE.xCon(i).c,2) == 1
                ind_tcon = [ind_tcon, i]
            end
        end
        Ic    = ind_tcon(spm_input('Which contrast?','!+1','m',{SwE.xCon(ind_tcon).name}));
        TITLE = {Cplot SwE.xCon(Ic).name};

        % select session and trial if
        %------------------------------------------------------------------
    case {'Event-related responses','Parametric responses','Volterra Kernels'}

        % get session
        %------------------------------------------------------------------
        s     = length(Sess);
        if  s > 1
            s = spm_input('which session','+1','n1',1,s);
        end

        % effect names
        %------------------------------------------------------------------
        switch Cplot
            case 'Volterra Kernels'
                u = length(Sess(s).Fc);
            otherwise
                u = length(Sess(s).U);
        end
        Uname = {};
        for i = 1:u
            Uname{i} = Sess(s).Fc(i).name;
        end

        % get effect
        %------------------------------------------------------------------
        str   = sprintf('which effect');
        u     = spm_input(str,'+1','m',Uname);

        % bin size
        %------------------------------------------------------------------
        dt    = SwE.xBF.dt;

end

spm('Pointer','Watch');

%-Extract filtered and whitened data from files
%==========================================================================
try
    y = spm_data_read(SwE.xY.VY, 'xyz', XYZ);
catch
    try
        % remap files in SPM.xY.P if SPM.xY.VY is no longer valid
        %------------------------------------------------------------------
        SwE.xY.VY = spm_data_hdr_read(SwE.xY.P);
        y = spm_data_read(SwE.xY.VY, 'xyz', XYZ);
        
    catch
        % data has been moved or renamed
        %------------------------------------------------------------------
        choice = questdlg({'Original data have been moved or renamed',...
            'How to proceed next?'},...
            [mfilename ': data files missing...'],...
            'Specify','Search','Ignore','Ignore');
        
        switch choice
            case 'Specify'
                [SwE.xY.P,sts] = ...
                    spm_select(numel(SwE.xY.VY),'image','Select images');
                if ~sts
                    [Y,y,beta,Bcov] = deal([]);
                    spm('Pointer','Arrow');
                    return;
                end
                SwE.xY.VY = spm_data_hdr_read(SwE.xY.P);
                for i = 1:numel(SwE.xY.VY)
                    SwE.xY.VY(i).pinfo(1:2,:) = ...
                        SwE.xY.VY(i).pinfo(1:2,:)*SwE.xGX.gSF(i);
                end
                y = spm_data_read(SwE.xY.VY, 'xyz', XYZ);
            case 'Search'
                SwE.xY.VY = spm_check_filename(SwE.xY.VY);
                y = spm_data_read(SwE.xY.VY, 'xyz', XYZ);
            otherwise
                y = [];
        end
    end
end

XYZstr = sprintf(' at [%g, %g, %g]',xyz);



%-Get parameter and hyperparameter estimates
%==========================================================================

%-Parameter estimates:   beta = xX.pKX*xX.K*y;
%----------------------------------------------------------------------
beta  = spm_data_read(SwE.Vbeta, 'xyz', XYZ);

%-Compute residuals
%--------------------------------------------------------------------------
if isempty(y)

    % make R = NaN so it will not be plotted
    %----------------------------------------------------------------------
    R   = NaN(size(SwE.xX.X,1),1);

else
    
    % residuals (non-whitened)
    %----------------------------------------------------------------------
    R   = y-SwE.xX.X*beta;

end   
it = 0;
it2 = 0;
Co = SwE.xCon(Ic).c;
Bcov = zeros(size(Co,2)*(size(Co,2)+1)/2,1);
% detect the indices of the betas of interest
if size(Co,2)==1
    ind = find(Co ~= 0);
else
    ind = find(any(Co'~=0));
end
for j = 1:size(Co,1)
    for jj = j:size(Co,1)
        it = it + 1;
        if any(j == ind) & any(jj == ind)
            it2 = it2+1;
            weight = Co(j,:)'*Co(jj,:);
            if (j~=jj)
                weight = weight + weight';
            end
            weight = weight(tril(ones(size(Co,2)))==1);
            Bcov = Bcov + weight * spm_data_read(SwE.Vcov_beta(it), 'xyz', XYZ);
        end
    end
end
CI    = spm_invTcdf(1-0.025,spm_data_read(SwE.xCon(Ic).Vedf, 'xyz', XYZ));

spm('Pointer','Arrow');

%-Plot
%==========================================================================

%-Colour specifications and index;
%--------------------------------------------------------------------------
Col   = [0 0 0; .8 .8 .8; 1 .5 .5];

switch Cplot

    %-Plot parameter estimates
    %======================================================================
    case 'Contrast estimates and 95% C.I.'

        % compute contrast of parameter estimates and 95% C.I.
        %------------------------------------------------------------------

        cbeta = Co'*beta;
        SE    = sqrt(Bcov);
        CI    = CI*SE;

        contrast.contrast      = cbeta;
        contrast.standarderror = SE;
        contrast.interval      = 2*CI;
        assignin('base','contrast',contrast)

        % bar chart
        %------------------------------------------------------------------
        figure(Fgraph)
        subplot(2,1,2)
        cla
        hold on

        % estimates
        %------------------------------------------------------------------
        h     = bar(cbeta);
        set(h,'FaceColor',Col(2,:))

        % standard error
        %------------------------------------------------------------------
        for j = 1:length(cbeta)
            line([j j],([CI(j) 0 - CI(j)] + cbeta(j)),...
                'LineWidth',6,'Color',Col(3,:))
        end

        title(TITLE,'FontSize',12)
        xlabel('contrast')
        ylabel(['contrast estimate',XYZstr])
        set(gca,'XLim',[0.4 (length(cbeta) + 0.6)])
        hold off

        % set Y to empty so outputs are assigned
        %------------------------------------------------------------------
        Y = [];

        
    %-All fitted effects or selected effects
    %======================================================================
    case 'Fitted responses'

        % predicted or adjusted response
        %------------------------------------------------------------------
        str   = 'predicted or adjusted response?';
        if spm_input(str,'!+1','b',{'predicted','adjusted'},[1 0]);

            % fitted (predicted) data (Y = X1*beta)
            %--------------------------------------------------------------
            Y = SwE.xX.X*SwE.xCon(Ic).c*pinv(SwE.xCon(Ic).c)*beta;
        else

            % fitted (corrected)  data (Y = X1o*beta)
            %--------------------------------------------------------------
            Y = spm_FcUtil('Yc',SwE.xCon(Ic),SwE.xX.xKXs,beta);

        end

        % adjusted data
        %------------------------------------------------------------------
        y     = Y + R;

        % get ordinates
        %------------------------------------------------------------------
        Xplot = {'an explanatory variable',...
                 'scan or time',...
                 'a user specified ordinate'};
        Cx    = spm_input('plot against','!+1','m',Xplot);

        % an explanatory variable
        %------------------------------------------------------------------
        if     Cx == 1

            str  = 'Which explanatory variable?';
            i    = spm_input(str,'!+1','m',SwE.xX.name);
            x    = SwE.xX.X(:,i);
            XLAB = SwE.xX.name{i};

        % scan or time
        %------------------------------------------------------------------
        elseif Cx == 2

            if isfield(SwE.xY,'RT')
                x    = SwE.xY.RT*[1:size(Y,1)]';
                XLAB = 'time {seconds}';
            else
                x    = [1:size(Y,1)]';
                XLAB = 'scan number';
            end

        % user specified
        %------------------------------------------------------------------
        elseif Cx == 3

            x    = spm_input('enter ordinate','!+1','e','',size(Y,1));
            XLAB = 'ordinate';

        end

        % plot
        %------------------------------------------------------------------
        figure(Fgraph)
        subplot(2,1,2)
        cla
        hold on
        [p q] = sort(x);
        if all(diff(x(q)))
            plot(x(q),Y(q),'LineWidth',4,'Color',Col(2,:));
            plot(x(q),y(q),':','Color',Col(1,:));
            plot(x(q),y(q),'.','MarkerSize',8, 'Color',Col(3,:));

        else
            plot(x(q),Y(q),'.','MarkerSize',16,'Color',Col(1,:));
            plot(x(q),y(q),'.','MarkerSize',8, 'Color',Col(2,:));
            xlim = get(gca,'XLim');
            xlim = [-1 1]*diff(xlim)/4 + xlim;
            set(gca,'XLim',xlim)

        end
        title(TITLE,'FontSize',12)
        xlabel(XLAB)
        ylabel(['response',XYZstr])
        legend('fitted','plus error')
        hold off
        
        
    %-Modeling evoked responses based on Sess
    %======================================================================
    case 'Event-related responses'

        % get plot type
        %------------------------------------------------------------------
        Rplot   = { 'fitted response and PSTH',...
            'fitted response and 95% C.I.',...
            'fitted response and adjusted data'};

        if isempty(y)
            TITLE = Rplot{2};
        else
            TITLE = Rplot{spm_input('plot in terms of','+1','m',Rplot)};
        end

        % plot
        %------------------------------------------------------------------
        switch TITLE
            case 'fitted response and PSTH'
                % build a simple FIR model subpartition (X); bin size = TR
                %----------------------------------------------------------
                BIN         = SwE.xY.RT;
                %BIN         = max(2,BIN);
                xBF         = SwE.xBF;
                U           = Sess(s).U(u);
                U.u         = U.u(:,1);
                xBF.name    = 'Finite Impulse Response';
                xBF.order   = round(32/BIN);
                xBF.length  = xBF.order*BIN;
                xBF         = spm_get_bf(xBF);
                BIN         = xBF.length/xBF.order;
                X           = spm_Volterra(U,xBF.bf,1);
                k           = SwE.nscan(s);
                X           = X([0:(k - 1)]*SwE.xBF.T + SwE.xBF.T0 + 32,:);

                % place X in SPM.xX.X
                %----------------------------------------------------------
                jX          = Sess(s).row;
                iX          = Sess(s).col(Sess(s).Fc(u).i);
                iX0         = [1:size(SwE.xX.X,2)];
                iX0(iX)     = [];
                X           = [X SwE.xX.X(jX,iX0)];
                X           = SwE.xX.W(jX,jX)*X;
                X           = [X SwE.xX.K(s).X0];

                % Re-estimate to get PSTH and CI
                %----------------------------------------------------------
                j           = xBF.order;
                xX          = spm_sp('Set',X);
                pX          = spm_sp('x-',xX);
                PSTH        = pX*y(jX);
                res         = spm_sp('r',xX,y(jX));
                df          = size(X,1) - size(X,2);
                bcov        = pX*pX'*sum(res.^2)/df;
                PSTH        = PSTH(1:j)/dt;
                PST         = [1:j]*BIN - BIN/2;
                PCI         = CI*sqrt(diag(bcov(1:j,(1:j))))/dt;
        end

        % basis functions and parameters
        %------------------------------------------------------------------
        X     = SwE.xBF.bf/dt;
        x     = ([1:size(X,1)] - 1)*dt;
        j     = Sess(s).col(Sess(s).Fc(u).i(1:size(X,2)));
        B     = beta(j);

        % fitted responses with standard error
        %------------------------------------------------------------------
        Y     = X*B;
        CI    = CI*sqrt(diag(X*Bcov(j,j)*X'));

        % peristimulus times and adjusted data (y = Y + R)
        %------------------------------------------------------------------
        pst   = Sess(s).U(u).pst;
        bin   = round(pst/dt);
        q     = find((bin >= 0) & (bin < size(X,1)));
        y     = R(Sess(s).row(:));
        pst   = pst(q);
        y     = y(q) + Y(bin(q) + 1);

        % plot
        %------------------------------------------------------------------
        figure(Fgraph)
        subplot(2,1,2)
        hold on
        switch TITLE

            case 'fitted response and PSTH'
                %----------------------------------------------------------
                errorbar(PST,PSTH,PCI)
                plot(PST,PSTH,'LineWidth',4,'Color',Col(2,:))
                plot(x,Y,'-.','Color',Col(3,:))

            case 'fitted response and 95% C.I.'
                %----------------------------------------------------------
                plot(x,Y,'Color',Col(2,:),'LineWidth',4)
                plot(x,Y + CI,'-.',x,Y - CI,'-.','Color',Col(1,:))

            case 'fitted response and adjusted data'
                %----------------------------------------------------------
                plot(x,Y,'Color',Col(2,:),'LineWidth',4)
                plot(pst,y,'.','Color',Col(3,:))

        end

        % label
        %------------------------------------------------------------------
        [i j] = max(Y);
        text(ceil(1.1*x(j)),i,Sess(s).Fc(u).name,'FontSize',8);
        title(TITLE,'FontSize',12)
        xlabel('peristimulus time {secs}')
        ylabel(['response',XYZstr])
        hold off


    %-Modeling evoked responses based on Sess
    %======================================================================
    case 'Parametric responses'

        % return gracefully if no parameters
        %------------------------------------------------------------------
        if ~Sess(s).U(u).P(1).h, return, end

        % basis functions
        %------------------------------------------------------------------
        bf    = SwE.xBF.bf;
        pst   = ([1:size(bf,1)] - 1)*dt;

        % orthogonalised expansion of parameteric variable
        %------------------------------------------------------------------
        str   = 'which parameter';
        p     = spm_input(str,'+1','m',{Sess(s).U(u).P.name});
        P     = Sess(s).U(u).P(p).P;
        q     = [];
        for i = 0:Sess(s).U(u).P(p).h;
            q = [q P.^i];
        end
        q     = spm_orth(q);

        % parameter estimates for this effect
        %------------------------------------------------------------------
        B     = beta(Sess(s).col(Sess(s).Fc(u).i));

        % reconstruct trial-specific responses
        %------------------------------------------------------------------
        Y     = zeros(size(bf,1),size(q,1));
        uj    = Sess(s).U(u).P(p).i;
        for i = 1:size(P,1)
            U      = sparse(1,uj,q(i,:),1,size(Sess(s).U(u).u,2));
            X      = kron(U,bf);
            Y(:,i) = X*B;
        end
        [P j] = sort(P);
        Y     = Y(:,j);

        % plot
        %------------------------------------------------------------------
        figure(Fgraph)
        subplot(2,2,3)
        surf(pst,P,Y')
        shading flat
        title(Sess(s).U(u).name{1},'FontSize',12)
        xlabel('PST {secs}')
        ylabel(Sess(s).U(u).P(p).name)
        zlabel(['responses',XYZstr])
        axis square

        % plot
        %------------------------------------------------------------------
        subplot(2,2,4)
        [j i] = max(mean(Y,2));
        plot(P,Y(i,:),'LineWidth',4,'Color',Col(2,:))
        str   = sprintf('response at %0.1fs',i*dt);
        title(str,'FontSize',12)
        xlabel(Sess(s).U(u).P(p).name)
        axis square
        grid on


    %-Modeling evoked responses based on Sess
    %======================================================================
    case 'Volterra Kernels'

        % Parameter estimates and basis functions
        %------------------------------------------------------------------
        bf    = SwE.xBF.bf/dt;
        pst   = ([1:size(bf,1)] - 1)*dt;

        % second order kernel
        %------------------------------------------------------------------
        if u > length(Sess(s).U)

            % Parameter estimates and kernel
            %--------------------------------------------------------------
            B     = beta(Sess(s).col(Sess(s).Fc(u).i));
            i     = 1;
            Y     = 0;
            for p = 1:size(bf,2)
                for q = 1:size(bf,2)
                    Y = Y + B(i)*bf(:,p)*bf(:,q)';
                    i = i + 1;
                end
            end

            % plot
            %--------------------------------------------------------------
            figure(Fgraph)
            subplot(2,2,3)
            imagesc(pst,pst,Y)
            axis xy
            axis image

            title('2nd order Kernel','FontSize',12);
            xlabel('peristimulus time {secs}')
            ylabel('peristimulus time {secs}')

            subplot(2,2,4)
            plot(pst,Y)
            axis square
            grid on

            title(Sess(s).Fc(u).name,'FontSize',12);
            xlabel('peristimulus time {secs}')


        % first  order kernel
        %------------------------------------------------------------------
        else
            B = beta(Sess(s).col(Sess(s).Fc(u).i(1:size(bf,2))));
            Y = bf*B;

            % plot
            %--------------------------------------------------------------
            figure(Fgraph)
            subplot(2,1,2)
            plot(pst,Y)
            grid on
            axis square

            title({'1st order Volterra Kernel' Sess(s).Fc(u).name},...
                'FontSize',12);
            xlabel('peristimulus time {secs}')
            ylabel(['impulse response',XYZstr])
        end

end


% Turn hold button off - this will alert the user to press it again
%--------------------------------------------------------------------------
try
    set(get(gcbo,'Userdata'),'Value',0);
catch
end


%-call Plot UI
%--------------------------------------------------------------------------
spm_results_ui('PlotUi',gca)
