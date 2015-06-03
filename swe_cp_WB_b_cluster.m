function swe_cp_WB_b_cluster(SwE,nBootstrap)
 
%-Say hello
%--------------------------------------------------------------------------
Finter = spm('CreateIntWin','off');
set(Finter,'name','Wild Bootstrap');
set(Finter,'vis','on')

%-Get SwE.mat[s] if necessary
%--------------------------------------------------------------------------
if nargin == 0
    P     = cellstr(spm_select(Inf,'^SwE\.mat$','Select SwE.mat[s]'));
    for i = 1:length(P)
        swd     = fileparts(P{i});
        load(fullfile(swd,'SwE.mat'));
        SwE.swd = swd;
        swe_cp_WB_b_cluster(SwE,nBootstrap);
    end
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
    SwE.WB.VYR;
catch
    spm('alert!','Please assign data to this design', mfilename);
    spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
    return
end

%==========================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%==========================================================================

%-MAXMEM is the maximum amount of data processed at a time (bytes)
%--------------------------------------------------------------------------
MAXMEM = spm_get_defaults('stats.maxmem');
mmv    = MAXMEM/8/SwE.nscan;
blksz  = ceil(mmv);                             %-block size
S      = size(SwE.xVol.XYZ,2); %- number of in-mask voxels
nbch   = ceil(S/ blksz);                                %-# blocks

 
%-Initialise XYZ matrix of in-mask voxel co-ordinates (real space)
%--------------------------------------------------------------------------
XYZ   = SwE.xVol.XYZ;
resamplingMatrix = NaN(SwE.nscan,1);
for iS = 1:SwE.Subj.nSubj
    resamplingMatrix(SwE.Subj.iSubj == SwE.Subj.uSubj(iS)) = binornd(1,0.5);
end
resamplingMatrix(resamplingMatrix == 0) = -1;

%-Cycle over bunches blocks within planes to avoid memory problems
%==========================================================================
str   = sprintf('Parameter estimation\nBootsrap # %i', nBootstrap);
spm_progress_bar('Init',100,str,'');
maxScore = nan;
minScore = nan;
conWB = SwE.WB.con;
nSizeCon = size(conWB,1);
rankCon = rank(conWB);

% only for the HCP data
edf = 79;

% score multiplicator value if contrast of rank > 1
if (nSizeCon > 1)
  scoreMult = (edf-rankCon+1)/edf/rankCon;
end

% activated voxels for cluster-wise inference
if (SwE.WB.clusterWise == 1)
  activatedVoxels = false(1,S);
  if (SwE.WB.stat == 'T')
    activatedVoxelsNeg = false(1,S);
  end
end


for bch = 1:nbch                     %-loop over blocks
    
  if bch ~= nbch
      index = (1+(bch-1)*blksz) : (bch * blksz);
      count = (bch-1)*blksz;
  else
      index = (1+(bch-1)*blksz) : size(XYZ,2);
      blksz = length(index);
  end
  %-Print progress information in command window
  %------------------------------------------------------------------
   str = sprintf('Block %3d/%-3d',bch,nbch);

  if bch == 1
      str2 = '';
  else
      str2 = repmat(sprintf('\b'),1,72);
  end
  fprintf('%s%-40s: %30s',str2,str,' ');


  %-Get data & construct analysis mask
  %=================================================================
  fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...read data')

  %-Get the data in mask, compute threshold & implicit masks
  %------------------------------------------------------------------
  Y_b = spm_get_data(SwE.WB.VYR, XYZ(:,index),false) + ...
      spm_get_data(SwE.WB.VResR, XYZ(:,index),false) .* repmat(resamplingMatrix,1,blksz);

  %==================================================================
  %-Proceed with General Linear Model (if there are voxels)
  %==================================================================

  %-General linear model: Ordinary least squares estimation
  %--------------------------------------------------------------
  fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...estimation');%-#

  beta  = SwE.WB.pX * Y_b;                     %-Parameter estimates
  if (SwE.WB.RSwE == 1)
    restmp = Y_b - SwE.WB.tmpR2 * beta;
  else
    restmp = Y_b - SwE.xX.X * beta;
  end
  clear Y_b
  if SwE.SS == 4 % SC2
    if (SwE.WB.RSwE == 1)
      resR = zeros(size(restmp));
      for i = 1:SwE.Subj.nSubj
          resR(SwE.Subj.iSubj==SwE.Subj.uSubj(i),:) = SwE.WB.corrR{i} *...
              restmp(SwE.Subj.iSubj==SwE.Subj.uSubj(i),:);
      end
    else% unrestricted residuals if needed
      resR = zeros(size(restmp));
      for i = 1:SwE.Subj.nSubj
          resR(SwE.Subj.iSubj==SwE.Subj.uSubj(i),:) = SwE.WB.corr{i} *...
              restmp(SwE.Subj.iSubj==SwE.Subj.uSubj(i),:);
      end
    end
  else
    if (SwE.WB.RSwE == 1)
      resR  = restmp.* repmat(SwE.WB.corrR, 1, blksz);
    else
      resR  = restmp.* repmat(SwE.WB.corr, 1, blksz);
    end
  end
  clear restmp

  %-Estimation of the data variance-covariance components (modified SwE) 
  %-SwE estimation (classic version)
  %--------------------------------------------------------------           
  if isfield(SwE.type,'modified')
      Cov_vis=zeros(SwE.Vis.nCov_vis,blksz);
      for i = SwE.WB.Ind_Cov_vis_diag
          Cov_vis(i,:) = mean(resR(SwE.WB.Flagk(i,:),:).^2);
      end
      for i = SwE.WB.Ind_Cov_vis_off_diag
          Cov_vis(i,:)= sum(resR(SwE.WB.Flagk(i,:),:).*resR(SwE.WB.Flagkk(i,:),:)).*...
              sqrt(Cov_vis(SwE.WB.Ind_Cov_vis_diag(SwE.WB.Ind_corr_diag(i,1)),:).*...
              Cov_vis(SwE.WB.Ind_Cov_vis_diag(SwE.WB.Ind_corr_diag(i,2)),:)./...
              sum(resR(SwE.WB.Flagk(i,:),:).^2)./...
              sum(resR(SwE.WB.Flagkk(i,:),:).^2));
      end
      %NaN may be produced in cov. estimation when one correspondant
      %variance are = 0, so set them to 0
      Cov_vis(isnan(Cov_vis))=0;
      %need to check if the eigenvalues of Cov_vis matrices are >=0
      for g = 1:SwE.Gr.nGr
          for iVox = 1:blksz
              tmp = zeros(SwE.Vis.nVis_g(g));
              tmp(tril(ones(SwE.Vis.nVis_g(g)))==1) = Cov_vis(SwE.WB.iGr_Cov_vis_g==g,iVox);
              tmp = tmp + tmp' - diag(diag(tmp));
              [V D] = eig(tmp);
              if any (diag(D)<0) %Bug corrected (BG - 19/09/13)
                  D(D<0) = 0;
                  tmp = V * D * V';
                  Cov_vis(SwE.WB.iGr_Cov_vis_g==g,iVox) = tmp(tril(ones(SwE.Vis.nVis_g(g)))==1); %Bug corrected (BG - 19/09/13)
              end
          end
      end       
  end

  if (nSizeCon == 1)
    if (SwE.WB.stat == 'T')
      score = (conWB * beta) ./ sqrt(SwE.WB.weightR * Cov_vis);
      minScore = min(minScore, min(score));
    else
      score = (conWB * beta).^2 ./ (SwE.WB.weightR * Cov_vis);
    end
    maxScore = max(maxScore, max(score));

    if (SwE.WB.clusterWise == 1)
      activatedVoxels(index) = score >= SwE.WB.primaryThreshold;
      if (SwE.WB.stat == 'T')
      	activatedVoxelsNeg(index) = score <= -SwE.WB.primaryThreshold;
      end
    end
  else
    % need to loop at every voxel
    cCovBc = SwE.WB.weightR * Cov_vis;
    cBeta = conWB * beta;
    score = zeros(1, blksz);
    for iVox = 1:blksz
      cCovBc_vox = zeros(nSizeCon);
      cCovBc_vox(tril(ones(nSizeCon))==1) = cCovBc(:,iVox);
      cCovBc_vox = cCovBc_vox + cCovBc_vox' - diag(diag(cCovBc_vox));
      score(iVox) = cBeta(:,iVox)' / cCovBc_vox * cBeta(:,iVox);
      score(iVox) = scoreMult * score(iVox); % need to multiply the score as described by the test
    end
    % save cluster information is needed
    if (SwE.WB.clusterWise == 1)
      activatedVoxels(index) = score >= SwE.WB.primaryThreshold;
    end
    maxScore = max(maxScore, max(score));
  end   
            
%-Report progress
%----------------------------------------------------------------------
fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done');   
spm_progress_bar('Set', 100 * bch / nbch);

end % (bch)

% compute the max cluster size if needed (so many ways this can be
% done... Not sure this solution is the best)
if (SwE.WB.clusterWise == 1)
	LocActivatedVoxels = SwE.xVol.XYZ(:,activatedVoxels);
	clusterAssignment = spm_clusters(LocActivatedVoxels);
  nCluster     = max(clusterAssignment);
  clusterSize = histc(clusterAssignment,1:nCluster);
  maxClusterSize = max(clusterSize);
  if (SwE.WB.stat == 'T')
    LocActivatedVoxelsNeg = SwE.xVol.XYZ(:,activatedVoxelsNeg);
    clusterAssignmentNeg = spm_clusters(LocActivatedVoxelsNeg);
    nClusterNeg     = max(clusterAssignmentNeg);
    clusterSizeNeg = histc(clusterAssignmentNeg,1:nClusterNeg);
    maxClusterSizeNeg = max(clusterSizeNeg);    
  end
end
fprintf('\n');                                                          %-#
spm_progress_bar('Clear')

%-Save remaining results files and analysis parameters
%==========================================================================
fprintf('%-40s: %30s','Saving results','...writing');


%-Save analysis Wmax and Wmin in W_b.mat file
%--------------------------------------------------------------------------
if (SwE.WB.clusterWise == 1)
  textMaxClusterSize =sprintf('maxClusterSize_WB_%05d', nBootstrap);
  if (SwE.WB.stat == 'T')
    textMaxClusterSizeNeg =sprintf('maxClusterSizeNeg_WB_%05d', nBootstrap);
  end
end
if (SwE.WB.stat == 'T')
  textMin = sprintf('minScore_WB_%05d', nBootstrap);
end
textMax = sprintf('maxScore_WB_%05d', nBootstrap);
if spm_matlab_version_chk('7') >=0
  if (SwE.WB.clusterWise == 1)
    save(textMaxClusterSize, 'maxClusterSize','-V6');
    if (SwE.WB.stat == 'T')
      save(textMaxClusterSizeNeg, 'maxClusterSizeNeg','-V6');
    end    
  end
  if (SwE.WB.stat == 'T')
    save(textMin,'minScore','-V6');
  end
  save(textMax,'maxScore','-V6');
else
  if (SwE.WB.clusterWise == 1)
    save(textMaxClusterSize, 'maxClusterSize');
    if (SwE.WB.stat == 'T')
      save(textMaxClusterSizeNeg, 'maxClusterSizeNeg');
    end
  end
  if (SwE.WB.stat == 'T')
    save(textMin,'minScore');
  end
  save(textMax,'maxScore');
end

%==========================================================================
%- E N D: Cleanup GUI
%==========================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#
%spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
fprintf('...use the results section for assessment\n\n')  
