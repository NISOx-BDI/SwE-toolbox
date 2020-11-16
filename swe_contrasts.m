function [SwE] = swe_contrasts(SwE,Ic)
% This function writes statistic & p value images for parametric analyses.
% =========================================================================
% Fills in SwE.xCon and writes the following images for a parametric
% analysis:
%   - swe_{vox|dpx|dat}_{T|F}stat_c{c#}.{nii|gii|mat}: T/F tatistic image
%   - swe_{vox|dpx|dat}_{zT|xF}stat_c{c#}.{nii|gii|mat}: Z/Chi square equivalent
%                                          statistic image
%   - swe_{vox|dpx|dat}_{T|F}stat_lp_c{c#}.{nii|gii|mat}: -log10 P image.
%   - swe_{vox|dpx|dat}_edf_c{c#}.{nii|gii|mat}: error degrees of freedom image.
% =========================================================================
% FORMAT [SwE] = SwE_contrasts(SwE,Ic)
% -------------------------------------------------------------------------
% Inputs/Outputs:
%   - SwE: SwE data structure
%   - Ic:  Indices of xCon to compute
% =========================================================================
% Modified version of spm_contrasts adapted for the SwE toolbox
% By Bryan Guillaume
% Version Info:  $Format:%ci$ $Format:%h$

% Temporary SwE variable to check for any changes to SwE. We want to avoid
% always having to save SwE.mat unless it has changed, because this is
% slow. A side benefit is one can look at results with just read
% privileges.
%--------------------------------------------------------------------------
tmpSwE = SwE;

%-Get and change to results directory
%--------------------------------------------------------------------------
try
    cd(SwE.swd);
end

% For Wild Bootstrap we already made contrast images, we just need to
% record them, if we haven't already.
%--------------------------------------------------------------------------
if isfield(SwE, 'WB')
    if ~isfield(SwE, 'xCon') || isempty(SwE.xCon)
        SwE = swe_contrasts_WB(SwE);
    end
    return
end

%-Get contrast definitions (if available)
%--------------------------------------------------------------------------
try
    xCon = SwE.xCon;
catch
    xCon = [];
end

%-Set all contrasts by default
%--------------------------------------------------------------------------
if nargin < 2
    Ic   = 1:length(xCon);
end

%-Check data format
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
            metadata = {};
        end
    else
        file_ext = spm_file_ext;
        file_data_type = 'vox';
        metadata = {};
    end
end

%-Map parameter files
%--------------------------------------------------------------------------

%-OLS estimators and covariance estimates
%--------------------------------------------------------------------------
Vbeta = SwE.Vbeta;
Vcov_beta = SwE.Vcov_beta;
dof_type = SwE.dof.dof_type;
if dof_type == 1
    Vcov_beta_g = SwE.Vcov_beta_g;
end
if dof_type>1
    Vcov_vis = SwE.Vcov_vis;
end

%-Compute & store contrast parameters, contrast/ESS images, & SwE images
%==========================================================================
spm('Pointer','Watch')
XYZ   = SwE.xVol.XYZ;
S=size(XYZ,2);
for i = 1:length(Ic)

    %-Canonicalise contrast structure with required fields
    %----------------------------------------------------------------------
    ic = Ic(i);
    %-Write contrast images?
    %======================================================================
    if isempty(xCon(ic).Vspm)
        if ~isMat
          Q = cumprod([1,SwE.xVol.DIM(1:2)'])*XYZ - ...
            sum(cumprod(SwE.xVol.DIM(1:2)'));
        end
        Co=xCon(ic).c;
        nBeta = size(Co,1);
        nSizeCon = size(Co,2);
        xCon(ic).eidf=rank(Co);
        % detect the indices of the betas of interest
        if nSizeCon==1
            ind = find(Co ~= 0);
        else
            ind = find(any(Co'~=0));
        end
        nCov_beta = (nBeta+1)*nBeta/2;

        % if ".mat" format, load data now
        if isMat
          beta = importdata(Vbeta);
          S = size(beta,2);
        end

        %-Compute contrast
        %------------------------------------------------------
        fprintf('\t%-32s: %30s',sprintf('contrast image %2d',ic),...
            '...computing');                                %-#
        str   = 'contrast computation';
        swe_progress_bar('Init',100,str,'');
        % if the Co is a vector, then create Co * Beta (Vcon)
        if nSizeCon==1
            if ~isMat
              V      = Vbeta(ind);
            end
            cBeta     = zeros(1,S);
            for j=1:numel(ind)
              if isMat
               	cBeta = cBeta + Co(ind(j)) * beta(ind(j),:);
              else
                cBeta = cBeta + Co(ind(j)) * swe_data_read(V(j),'xyz',XYZ);
              end
                swe_progress_bar('Set',100*(j/numel(ind)));
            end
            swe_progress_bar('Clear')

            if isMat
              %-save contrasted beta
              %------------------------------------------------------
              xCon(ic).Vcon = sprintf('swe_%s_con_c%02d%s',file_data_type,ic,file_ext);
              save(xCon(ic).Vcon, 'cBeta')
            else
              %-Prepare handle for contrast image
              %------------------------------------------------------
              xCon(ic).Vcon = swe_data_hdr_write(sprintf('swe_%s_con_c%02d%s',file_data_type,ic,file_ext),...
                                                  SwE.xVol.DIM', SwE.xVol.M,...
                                                  sprintf('SwE contrast - %d: %s',ic,xCon(ic).name),...
                                                  metadata);
              %-Write image
              %------------------------------------------------------
              tmp = zeros(SwE.xVol.DIM');
              tmp(Q) = cBeta;
              xCon(ic).Vcon = swe_data_write(xCon(ic).Vcon, tmp);

              clear tmp
            end
            if ~isMat
              fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
                '...written %s', xCon(ic).Vcon.fname))
            else
              fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
                '...written %s', xCon(ic).Vcon))
            end
        else
            if ~isMat
              V      = Vbeta(ind);
            end
            cBeta     = zeros(nSizeCon,S);
            for j=1:numel(ind)
              if isMat
               	cBeta = cBeta + Co(ind(j),:)' * beta(ind(j),:);
              else
                cBeta = cBeta + Co(ind(j),:)' * swe_data_read(V(j),'xyz',XYZ);
              end
              swe_progress_bar('Set',100*(j/numel(ind)));
            end
            swe_progress_bar('Clear')
            if ~isMat
              fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
                '...computed'))
            else
              fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
                '...computed'))
            end
        end

        % clear beta for memory
        if isMat
          clear beta
        end
        %-Write inference SwE
        %======================================================================

        %-compute the contrasted beta covariances and edof for the contrast
        switch(xCon(ic).STAT)
          case 'T'
            eSTAT = 'Z';
          case 'F'
            eSTAT = 'X';
        end
        fprintf('\t%-32s: %30s',sprintf('spm{%c} image %2d',eSTAT,ic),...
            '...computing');                                %-#
        str   = 'contrasted beta covariance computation';
        swe_progress_bar('Init',100,str,'');

        it = 0;
        it2 = 0;
        cCovBc = zeros(nSizeCon*(nSizeCon+1)/2,S);
        if dof_type == 1
            cCovBc_g = zeros(nSizeCon*(nSizeCon+1)/2,S,SwE.Gr.nGr);
        else
            indSubDesignMatrices = SwE.dof.iBeta_dof(ind);
            subjectsInvolved = [];
            for iIndSubDesignMatrices = indSubDesignMatrices
              subjectsInvolved = [subjectsInvolved; SwE.Subj.iSubj(SwE.dof.iGr_dof == iIndSubDesignMatrices)];
            end
            subjectsInvolved = unique(subjectsInvolved);
            indSubjInvolved = nan(length(subjectsInvolved),1);
            % convert into in
            for iSubjInvolved = 1:length(subjectsInvolved)
              indSubjInvolved(iSubjInvolved) = find(SwE.Subj.uSubj == subjectsInvolved(iSubjInvolved));
            end
            xCon(ic).edf = sum(SwE.dof.edof_Subj(indSubjInvolved));
        end

        % load .mat file(s) if this is the format
        if isMat
          cov_beta = importdata(Vcov_beta);
          if dof_type == 1
            cov_beta_g = importdata(Vcov_beta_g);
          end
        end

        for j = 1:nBeta
            for jj = j:nBeta
                it = it + 1;
                if any(j == ind) && any(jj == ind)
                    it2 = it2+1;
                    weight = Co(j,:)'*Co(jj,:);
                    if (j~=jj) %was wrong (BG - 13/09/13)
                        weight = weight + weight';
                    end
                    weight = weight(tril(ones(nSizeCon))==1);
                    if isMat
                      cCovBc = cCovBc + weight * cov_beta(it,:);
                    else
                      cCovBc = cCovBc + weight * swe_data_read(Vcov_beta(it),'xyz',XYZ);
                    end
                    if dof_type == 1
                      for g = 1:SwE.Gr.nGr
                        if isMat
                          cCovBc_g(:,:,g) = cCovBc_g(:,:,g) + weight *...
                            reshape(cov_beta_g(g,it,:), 1, S);
                        else
                          cCovBc_g(:,:,g) = cCovBc_g(:,:,g) + weight *...
                            swe_data_read(Vcov_beta_g((g-1)*nCov_beta+it),'xyz',XYZ);
                        end
                        swe_progress_bar('Set',100*((it2-1+g/SwE.Gr.nGr)/length(ind)/(length(ind)+1)*2));
                      end
                    end
                    swe_progress_bar('Set',100*(it2/length(ind)/(length(ind)+1)*2));
                end
            end
        end
        swe_progress_bar('Clear')

        str   = 'spm computation';
        swe_progress_bar('Init',100,str,'');
        equivalentScore = zeros(1,S);
        % add output of uncorrected p-values
        uncP            = zeros(1,S);
        switch(xCon(ic).STAT)
            case 'T'                                 %-Compute spm{t} image
                %----------------------------------------------------------
                score = cBeta ./ sqrt(cCovBc);
                swe_progress_bar('Set',100*(0.1));
                switch dof_type
                    case 1
                        tmp = 0;
                        for g = 1:SwE.Gr.nGr
                            tmp = tmp + cCovBc_g(:,:,g).^2/SwE.dof.edof_Gr(g);
                            swe_progress_bar('Set',100*(g/SwE.Gr.nGr/10+0.1));
                        end
                        clear cCovBc_g
                        edf = cCovBc.^2 ./ tmp;
                        swe_progress_bar('Set',100*(0.2));
                        % transform into Z-scores image
                        if any(score>0) % avoid to run the following line when all Z are < 0 (BG - 22/08/2016)
                          uncP(score>0)             = spm_Tcdf(-score(score>0),edf(score>0));
                          equivalentScore(score>0)  = -swe_invNcdf(uncP(score>0));
                        end
                        if any(score<0) % avoid to run the following line when all Z are > 0(BG - 22/08/2016)
                          uncP(score<0)             = spm_Tcdf(score(score<0),edf(score<0));
                          equivalentScore(score<0)  = swe_invNcdf(uncP(score<0));
                          uncP(score<0)             = 1 - uncP(score<0);
                        end
                        %Z = -log10(1-spm_Tcdf(Z,edf)); %transfo into -log10(p)
                        swe_progress_bar('Set',100);
                    case 0
                        % transform into Z-scores image
                        if any(score>0) % avoid to run the following line when all Z are < 0 (BG - 22/08/2016)
                          uncP(score>0)             = spm_Tcdf(-score(score>0),xCon(ic).edf);
                          equivalentScore(score>0)  = -swe_invNcdf(uncP(score>0));
                        end
                        if any(score<0) % avoid to run the following line when all Z are > 0(BG - 22/08/2016)
                          uncP(score<0)             = spm_Tcdf(score(score<0),xCon(ic).edf);
                          equivalentScore(score<0)  = swe_invNcdf(uncP(score<0));
                          uncP(score<0)             = 1 - uncP(score<0);
                        end
                        % transform into -log10(p-values) image
                        %Z = -log10(1-spm_Tcdf(Z,xCon(ic).edf));
                        swe_progress_bar('Set',100);
                    case 2
                        CovcCovBc = 0;
                        if isMat
                          cov_vis = importdata(Vcov_vis);
                        end
                        for g = 1:SwE.Gr.nGr
                            Wg = kron(Co,Co)' * swe_duplication_matrix(nBeta) * SwE.Vis.weight(:,SwE.Vis.iGr_Cov_vis_g==g);
                            Wg = kron(Wg,Wg) * swe_duplication_matrix(SwE.Vis.nCov_vis_g(g));
                            if isMat
                            	CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(cov_vis(SwE.Vis.iGr_Cov_vis_g==g,:),SwE.dof.dofMat{g},1);
                            else
                              CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(swe_data_read(Vcov_vis(SwE.Vis.iGr_Cov_vis_g==g),'xyz',XYZ),SwE.dof.dofMat{g},1);
                            end
                            swe_progress_bar('Set',100*(0.1) + g*80/SwE.Gr.nGr);
                        end
                        clear Wg cov_vis
                        edf = 2 * cCovBc.^2 ./ CovcCovBc - 2;
                        clear CovcCovBc
                        if any(score>0) % avoid to run the following line when all Z are < 0 (BG - 22/08/2016)
                          uncP(score>0)             = spm_Tcdf(-score(score>0),edf(score>0));
                          equivalentScore(score>0)  = -swe_invNcdf(uncP(score>0));
                        end
                        if any(score<0) % avoid to run the following line when all Z are > 0(BG - 22/08/2016)
                          uncP(score<0)             = spm_Tcdf(score(score<0),edf(score<0));
                          equivalentScore(score<0)  = swe_invNcdf(uncP(score<0));
                          uncP(score<0)             = 1 - uncP(score<0);
                        end
                        %Z = -log10(1-spm_Tcdf(Z,edf)); %transfo into -log10(p)
                        swe_progress_bar('Set',100);
                    case 3
                        CovcCovBc = 0;
                        if isMat
                          cov_vis = importdata(Vcov_vis);
                        end
                        for g = 1:SwE.Gr.nGr
                          Wg = kron(Co,Co)' * swe_duplication_matrix(nBeta) * SwE.Vis.weight(:,SwE.Vis.iGr_Cov_vis_g==g);
                          Wg = kron(Wg,Wg) * swe_duplication_matrix(SwE.Vis.nCov_vis_g(g));
                          if isMat
                            CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(cov_vis(SwE.Vis.iGr_Cov_vis_g==g,:),SwE.dof.dofMat{g},2);
                          else
                            CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(swe_data_read(Vcov_vis(SwE.Vis.iGr_Cov_vis_g==g),'xyz',XYZ),SwE.dof.dofMat{g},2);
                          end
                          swe_progress_bar('Set',100*(0.1) + g*80/SwE.Gr.nGr);
                        end
                        clear Wg cov_vis
                        edf = 2 * cCovBc.^2 ./ CovcCovBc;
                        clear CovcCovBc
                        % transform into Z-scores image
                        if any(score>0) % avoid to run the following line when all Z are < 0 (BG - 22/08/2016)
                          uncP(score>0)             = spm_Tcdf(-score(score>0),edf(score>0));
                          equivalentScore(score>0)  = -swe_invNcdf(uncP(score>0));
                        end
                        if any(score<0) % avoid to run the following line when all Z are > 0(BG - 22/08/2016)
                          uncP(score<0)             = spm_Tcdf(score(score<0),edf(score<0));
                          equivalentScore(score<0)  = swe_invNcdf(uncP(score<0));
                          uncP(score<0)             = 1 - uncP(score<0);
                        end
                        %Z = -log10(1-spm_Tcdf(Z,edf)); %transfo into -log10(p)
                        swe_progress_bar('Set',100);
                end

            case 'F'                                 %-Compute spm{F} image
                %---------------------------------------------------------
                if nSizeCon==1
                    score = abs(cBeta ./ sqrt(cCovBc));
                    indNotNan = ~isnan(score);
                    swe_progress_bar('Set',100*(0.1));
                    switch dof_type
                        case 1
                            tmp = 0;
                            for g = 1:SwE.Gr.nGr
                                tmp = tmp + cCovBc_g(:,:,g).^2/SwE.dof.edof_Gr(g);
                                swe_progress_bar('Set',100*(g/SwE.Gr.nGr/10+0.1));
                            end
                            clear cCovBc_g
                            edf = cCovBc.^2 ./ tmp;
                            swe_progress_bar('Set',100*(3/4));
                            % transform into X-scores image
                            uncP(indNotNan) = spm_Tcdf(-abs(score(indNotNan)), edf(indNotNan));
                            equivalentScore(indNotNan) = (swe_invNcdf(uncP(indNotNan))).^2;
                            uncP(indNotNan) = 2 *  uncP(indNotNan);
                            % transform into -log10(p-values) image
                            %Z = -log10(1-spm_Fcdf(Z,1,edf));
                            swe_progress_bar('Set',100);
                        case 0
                            % transform into X-scores image
                            uncP(indNotNan) = spm_Tcdf(-abs(score(indNotNan)),xCon(ic).edf);
                            equivalentScore(indNotNan) = (swe_invNcdf(uncP(indNotNan))).^2;
                            uncP(indNotNan) = 2 *  uncP(indNotNan);
                           % transform into -log10(p-values) image
                            %Z = -log10(1-spm_Fcdf(Z,1, xCon(ic).edf));
                            swe_progress_bar('Set',100);
                        case 2
                            CovcCovBc = 0;
                            if isMat
                              cov_vis = importdata(Vcov_vis);
                            end
                            for g = 1:SwE.Gr.nGr
                                Wg = kron(Co,Co)' * swe_duplication_matrix(nBeta) * SwE.Vis.weight(:,SwE.Vis.iGr_Cov_vis_g==g);
                                Wg = kron(Wg,Wg) * swe_duplication_matrix(SwE.Vis.nCov_vis_g(g));
                                if isMat
                                  CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(cov_vis(SwE.Vis.iGr_Cov_vis_g==g,:),SwE.dof.dofMat{g},1);
                                else
                                  CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(swe_data_read(Vcov_vis(SwE.Vis.iGr_Cov_vis_g==g),'xyz',XYZ),SwE.dof.dofMat{g},1);
                                end
                                swe_progress_bar('Set',100*(g/SwE.Gr.nGr/10+0.1));
                            end
                            clear Wg cov_vis
                            edf = 2 * cCovBc.^2 ./ CovcCovBc - 2;
                            clear CovcCovBc
                            swe_progress_bar('Set',100*(3/4));
                            % transform into X-scores image
                            uncP(indNotNan) = spm_Tcdf(-abs(score(indNotNan)), edf(indNotNan));
                            equivalentScore(indNotNan) = (swe_invNcdf(uncP(indNotNan))).^2;
                            uncP(indNotNan) = 2 *  uncP(indNotNan);
                            % transform into -log10(p-values) image
                            %Z = -log10(1-spm_Fcdf(Z,1,edf));
                            swe_progress_bar('Set',100);
                        case 3
                            CovcCovBc = 0;
                            if isMat
                              cov_vis = importdata(Vcov_vis);
                            end
                            for g = 1:SwE.Gr.nGr
                                Wg = kron(Co,Co)' * swe_duplication_matrix(nBeta) * SwE.Vis.weight(:,SwE.Vis.iGr_Cov_vis_g==g);
                                Wg = kron(Wg,Wg) * swe_duplication_matrix(SwE.Vis.nCov_vis_g(g));
                                if isMat
                                  CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(cov_vis(SwE.Vis.iGr_Cov_vis_g==g,:),SwE.dof.dofMat{g},2);
                                else
                                  CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(swe_data_read(Vcov_vis(SwE.Vis.iGr_Cov_vis_g==g),'xyz',XYZ),SwE.dof.dofMat{g},2);
                                end
                                swe_progress_bar('Set',100*(g/SwE.Gr.nGr/10+0.1));
                            end
                            clear Wg cov_vis
                            edf = 2 * cCovBc.^2 ./ CovcCovBc;
                            % transform into X-scores image
                            uncP(indNotNan) = spm_Tcdf(-abs(score(indNotNan)), edf(indNotNan));
                            equivalentScore(indNotNan) = (swe_invNcdf(uncP(indNotNan))).^2;
                            uncP(indNotNan) = 2 *  uncP(indNotNan);
                            % transform into -log10(p-values) image
                            %Z = -log10(1-spm_Fcdf(Z,1,edf));
                            swe_progress_bar('Set',100);
                            clear CovcCovBc
                    end
                    % need to transform in F-score, not in absolute t-score
                    % corrected on 12/05/15 by BG
                    score = score.^2;
                else
                    score   = zeros(1,S);
                    if dof_type ~= 0
                        edf = zeros(1,S);
                    end
                    if dof_type == 2
                        CovcCovBc = 0;
                        if isMat
                          cov_vis = importdata(Vcov_vis);
                        end
                        for g = 1:SwE.Gr.nGr
                             Wg = kron(Co,Co)' * swe_duplication_matrix(nBeta) * SwE.Vis.weight(:,SwE.Vis.iGr_Cov_vis_g==g);
                             Wg = sum(kron(Wg,Wg)) * swe_duplication_matrix(SwE.Vis.nCov_vis_g(g));
                             if isMat
                               CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(cov_vis(SwE.Vis.iGr_Cov_vis_g==g,:),SwE.dof.dofMat{g},1);
                             else
                               CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(swe_data_read(Vcov_vis(SwE.Vis.iGr_Cov_vis_g==g),'xyz',XYZ),SwE.dof.dofMat{g},1);
                             end
                        end
                        clear cov_vis
                        edf = 2 * (sum(swe_duplication_matrix(nSizeCon)) * cCovBc).^2 ./ CovcCovBc - 2;
                    end
                    if dof_type == 3
                      CovcCovBc = 0;
                      if isMat
                        cov_vis = importdata(Vcov_vis);
                      end
                      tmp = eye(nSizeCon^2);
                      for g = 1:SwE.Gr.nGr
                        Wg = kron(Co,Co)' * swe_duplication_matrix(nBeta) * SwE.Vis.weight(:,SwE.Vis.iGr_Cov_vis_g==g);
                        % tmp is used to sum only the diagonal element
                        % this is useful to compute the trace as
                        % tr(A) = vec(I)' * vec(A)
                        Wg = tmp(:)' * (kron(Wg,Wg)) * swe_duplication_matrix(SwE.Vis.nCov_vis_g(g));
                        if isMat
                          CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(cov_vis(SwE.Vis.iGr_Cov_vis_g==g,:),SwE.dof.dofMat{g},2);
                        else
                          CovcCovBc = CovcCovBc + Wg * swe_vechCovVechV(swe_data_read(Vcov_vis(SwE.Vis.iGr_Cov_vis_g==g),'xyz',XYZ),SwE.dof.dofMat{g},2);
                        end
                      end
                      clear cov_vis
                      % note that tr(A^2) = vec(A)' * vec(A)
                      tmp = eye(nSizeCon);
                      edf = (sum(swe_duplication_matrix(nSizeCon)) * cCovBc.^2 +...
                        (tmp(:)' * swe_duplication_matrix(nSizeCon) * cCovBc).^2) ./ CovcCovBc;
                    end
                    % define a parameter to tell when to update progress
                    % bar only 80 times
                    updateEvery = round(S/80);
                    indNotNan = ~isnan(cCovBc(1,:));
                    for iVox=1:S
                      if indNotNan(iVox)
                        cCovBc_vox = zeros(nSizeCon);
                        cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
                        cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
                        score(iVox) = cBeta(:,iVox)' / cCovBc_vox * cBeta(:,iVox);
                        if (dof_type == 1)
                          tmp = 0;
                          for g = 1:SwE.Gr.nGr
                            cCovBc_g_vox = zeros(nSizeCon);
                            cCovBc_g_vox(tril(ones(nSizeCon))==1) = cCovBc_g(:,iVox,g);
                            cCovBc_g_vox = cCovBc_g_vox + cCovBc_g_vox' - diag(diag(cCovBc_g_vox));
                            tmp = tmp + (trace(cCovBc_g_vox^2) + (trace(cCovBc_g_vox))^2)/...
                              SwE.dof.edof_Gr(g);
                          end
                          edf(iVox)=(trace(cCovBc_vox^2) + (trace(cCovBc_vox))^2) / tmp;
                        end
                        % update progress_bar only approx 80 times
                        if (mod(iVox,updateEvery) == 0)
                          swe_progress_bar('Set',10 + 80 * (iVox/S));
                        end
                      end
                    end
                    if dof_type ~= 0
                        clear cCovBc_g
                        score = score .*(edf-xCon(ic).eidf+1)./edf/xCon(ic).eidf;
                        score(score < 0) = 0; % force negatif F-score to 0 (can happen for very low edf)
                        % transform into X-scores image
                        % Z2 = chi2inv(spm_Fcdf(Z,xCon(ic).eidf,edf-xCon(ic).eidf+1),1);
                        try % check if the user do not have the fcdf function or one with 'upper' option
                          uncP(score>1)             = fcdf(score(score>1),xCon(ic).eidf,edf(score>1)-xCon(ic).eidf+1,'upper'); % more accurate to look this way for high scores
                          equivalentScore(score>1)  = swe_invNcdf(uncP(score>1)/2).^2;
                          uncP(score<=1 & score > 0)            = 1 - fcdf(score(score<=1 & score > 0),xCon(ic).eidf,edf(score<=1 & score > 0)-xCon(ic).eidf+1);
                          equivalentScore(score<=1 & score > 0) = swe_invNcdf(uncP(score<=1 & score > 0)/2).^2;
                        catch
                          uncP(score>0)            = betainc((edf(score>0) - xCon(ic).eidf + 1)./(edf(score>0) - xCon(ic).eidf + 1 + xCon(ic).eidf * score(score>0)),(edf(score>0)-xCon(ic).eidf+1)/2, xCon(ic).eidf/2); % more accurate to look this way for high scores
                          equivalentScore(score>0) = swe_invNcdf(uncP(score>0)/2).^2;
%                             Z2 = swe_invNcdf(0.5 - spm_Fcdf(Z,xCon(ic).eidf, edf-xCon(ic).eidf+1)/2).^2;
                        end
                        uncP(score == 0)            = 0;
                        equivalentScore(score == 0) = 0;
                        % transform into -log10(p-values) image
                        %Z = -log10(1-spm_Fcdf(Z,xCon(ic).eidf,edf));
                    else
                        score = score *(xCon(ic).edf -xCon(ic).eidf+1)/xCon(ic).edf/xCon(ic).eidf;
                        score(score < 0) = 0; % force negatif F-score to 0 (can happen for very low edf)
                        % transform into X-scores image
                        %Z2 = chi2inv(spm_Fcdf(Z,xCon(ic).eidf,xCon(ic).edf-xCon(ic).eidf+1),1);
                        try % check if the user do not have the fcdf function or one with 'upper' options
                          uncP(score>1)             = fcdf(score(score>1),xCon(ic).eidf,xCon(ic).edf-xCon(ic).eidf+1,'upper');% more accurate to look this way for high score
                          equivalentScore(score>1)  = swe_invNcdf(uncP(score>1)/2).^2;
                          uncP(score<=1 & score > 0)            = 1 - fcdf(score(score<=1 & score > 0),xCon(ic).eidf,xCon(ic).edf-xCon(ic).eidf+1);
                          equivalentScore(score<=1 & score > 0) = swe_invNcdf(uncP(score<=1 & score > 0)/2).^2;
                        catch
                          uncP(score>0)            = betainc((xCon(ic).edf - xCon(ic).eidf + 1)./(xCon(ic).edf - xCon(ic).eidf + 1 + xCon(ic).eidf * score(score>0)),(xCon(ic).edf-xCon(ic).eidf+1)/2, xCon(ic).eidf/2);
                          equivalentScore(score>0) = swe_invNcdf(uncP(score>0)/2).^2;
%                           Z2(Z>0) = swe_invNcdf(0.5 - spm_Fcdf(Z,xCon(ic).eidf,xCon(ic).edf-xCon(ic).eidf+1)/2).^2;
                        end
                        uncP(score == 0)            = 0;
                        equivalentScore(score == 0) = 0;
                        % transform into -log10(p-values) image
                        %Z = -log10(1-spm_Fcdf(Z,xCon(ic).eidf,xCon(ic).edf));
                    end
                    swe_progress_bar('Set',100);
                end
        end
        luncP = -log10(uncP);
        swe_progress_bar('Clear')
        clear cCovBc cB tmp


        %-Write SwE - statistic images & edf image if needed
        %------------------------------------------------------------------
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...writing');      %-#

        equivalentScore (equivalentScore > realmax('single')) = realmax('single');
        equivalentScore (equivalentScore < -realmax('single')) = -realmax('single');

        if isMat
          xCon(ic).Vspm = sprintf('swe_%s_%c%cstat_c%02d%s',file_data_type,lower(eSTAT),xCon(ic).STAT,ic, file_ext);
          save(xCon(ic).Vspm, 'equivalentScore');
        else
          xCon(ic).Vspm = swe_data_hdr_write(sprintf('swe_%s_%c%cstat_c%02d%s',file_data_type,lower(eSTAT),xCon(ic).STAT,ic,file_ext),...
                                              SwE.xVol.DIM', SwE.xVol.M,...
                                              sprintf('spm{%c} - contrast %d: %s',eSTAT,ic,xCon(ic).name),...
                                              metadata);

          tmp           = zeros(SwE.xVol.DIM');
          tmp(Q)        = equivalentScore;
          xCon(ic).Vspm = swe_data_write(xCon(ic).Vspm,tmp);
        end
        clear tmp equivalentScore
        if isMat
          fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(ic).Vspm,'t')));
        else
          fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(ic).Vspm.fname,'t')));
        end

        fprintf('\t%-32s: %30s',sprintf('spm{%c} image %2d',xCon(ic).STAT,ic),...
        '...writing');

        if isMat
          xCon(ic).Vspm2 = sprintf('swe_%s_%cstat_c%02d%s',file_data_type,xCon(ic).STAT,ic,file_ext);
          save(xCon(ic).Vspm2, 'score');
        else
          xCon(ic).Vspm2 = struct(...
            'fname',  sprintf('swe_%s_%cstat_c%02d%s',file_data_type,xCon(ic).STAT,ic,file_ext),...
            'dim',    SwE.xVol.DIM',...
            'dt',     [spm_type('float32'), spm_platform('bigend')],...
            'mat',    SwE.xVol.M,...
            'pinfo',  [1,0,0]',...
            'descrip',sprintf('spm{%c} - contrast %d: %s',xCon(ic).STAT,ic,xCon(ic).name),...
            metadata{:});

          xCon(ic).Vspm2 = swe_data_hdr_write(sprintf('swe_%s_%cstat_c%02d%s',file_data_type,xCon(ic).STAT,ic,file_ext),...
                                              SwE.xVol.DIM', SwE.xVol.M,...
                                              sprintf('spm{%c} - contrast %d: %s',xCon(ic).STAT,ic,xCon(ic).name),...
                                              metadata);

          tmp           = zeros(SwE.xVol.DIM');
          tmp(Q)        = score;
          xCon(ic).Vspm2 = swe_data_write(xCon(ic).Vspm2,tmp);
        end
        clear tmp score
        if isMat
          fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(ic).Vspm2,'t')));   %-#
        else
           fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(ic).Vspm2.fname,'t')));   %-#
        end

        % save raw uncorrected p-values (new on 05/11/2017)
        fprintf('\t%-32s: %30s',sprintf('-log10(uncP) image %2d',ic),...
        '...writing');
        if isMat
          xCon(ic).VspmUncP = sprintf('swe_%s_%cstat_lp_c%02d%s',file_data_type,xCon(ic).STAT,ic,file_ext);
          save(xCon(ic).VspmUncP, 'luncP');
        else
          xCon(ic).VspmUncP = struct(...
            'fname',  sprintf('swe_%s_%cstat_lp_c%02d%s',file_data_type,xCon(ic).STAT,ic,file_ext),...
            'dim',    SwE.xVol.DIM',...
            'dt',     [spm_type('float32'), spm_platform('bigend')],...
            'mat',    SwE.xVol.M,...
            'pinfo',  [1,0,0]',...
            'descrip',sprintf('spm{%c} - contrast %d: %s','luncP',ic,xCon(ic).name),...
            metadata{:});

          xCon(ic).VspmUncP = swe_data_hdr_write(sprintf('swe_%s_%cstat_lp_c%02d%s',file_data_type,xCon(ic).STAT,ic,file_ext),...
                                                  SwE.xVol.DIM', SwE.xVol.M,...
                                                  sprintf('spm{%c} - contrast %d: %s','luncP',ic,xCon(ic).name),...
                                                  metadata);

          tmp           = zeros(SwE.xVol.DIM');
          tmp(Q)        = luncP;
          xCon(ic).VspmUncP = swe_data_write(xCon(ic).VspmUncP,tmp);
        end
        clear tmp uncP luncP
        if isMat
          fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(ic).VspmUncP,'t')));
        else
          fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(ic).VspmUncP.fname,'t')));
        end

        fprintf('\t%-32s: %30s',sprintf('edf image %2d',ic),...
        '...writing');

        if dof_type
          if isMat
            xCon(ic).Vedf = sprintf('swe_%s_edf_c%02d%s',file_data_type,ic,file_ext);
            save(xCon(ic).Vedf, 'edf');
          else
            xCon(ic).Vedf = struct(...
              'fname',  sprintf('swe_%s_edf_c%02d%s',file_data_type,ic,file_ext),...
              'dim',    SwE.xVol.DIM',...
              'dt',     [spm_type('float32'), spm_platform('bigend')],...
              'mat',    SwE.xVol.M,...
              'pinfo',  [1,0,0]',...
              'descrip',sprintf('SwE effective degrees of freedom - %d: %s',ic,xCon(ic).name),...
              metadata{:});

            xCon(ic).Vedf = swe_data_hdr_write(sprintf('swe_%s_edf_c%02d%s',file_data_type,ic,file_ext),...
                                                SwE.xVol.DIM', SwE.xVol.M,...
                                                sprintf('SwE effective degrees of freedom - %d: %s',ic,xCon(ic).name),...
                                                metadata);

            tmp = zeros(SwE.xVol.DIM');
            tmp(Q) = edf;
            xCon(ic).Vedf = swe_data_write(xCon(ic).Vedf,tmp);
          end
          clear tmp edf
          if isMat
            fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
              '...written %s',spm_str_manip(xCon(ic).Vedf,'t')))%-#
          else
            fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
              '...written %s',spm_str_manip(xCon(ic).Vedf.fname,'t')))%-#
          end
        end

    end % if isempty(xCon(ic).Vspm)

end % (for i = 1:length(Ic))
spm('Pointer','Arrow')

% place xCon back in SwE
%--------------------------------------------------------------------------
SwE.xCon = xCon;

% Check if SwE has changed. Save only if it has.
%--------------------------------------------------------------------------
if isOctave
    save('SwE.mat','SwE');
elseif spm_check_version('matlab','7') >=0
    save('SwE','SwE','-V6');
else
    save('SwE','SwE');
end
