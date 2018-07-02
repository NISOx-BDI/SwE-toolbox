% This function creates the missing 'xCon' needed for wild bootstrap SwE
% objects. Note: Whilst contrasts are created in 'swe_contrasts.m', in 
% 'swe_WB_contrasts.m' they are only recorded. They have already been
% created by 'swe_cp_WB.m'.
%
% Author: Tom Maullin (08/06/2018)
% =========================================================================

function [SwE] = swe_contrasts_WB(SwE)

    %-Get and change to results directory
    %----------------------------------------------------------------------
    try
        cd(SwE.swd);
    end
    
    %-Get file extension
    %----------------------------------------------------------------------
    [~,~,file_ext] = fileparts(SwE.xY.P{1});
    isMat          = strcmpi(file_ext,'.mat');

    if ~isMat
        isMeshData = spm_mesh_detect(SwE.xY.VY);
        if isMeshData
            file_ext = '.gii';
        else
            file_ext = spm_file_ext;
        end
    end
    
    %-Retrieve contrast details.
    %----------------------------------------------------------------------
    STAT = SwE.WB.stat;
    c = SwE.WB.con;
    
    %-Retrieve equivalent statistic details.
    %----------------------------------------------------------------------
    if strcmpi(STAT, 't')
        eSTAT = 'z';
    else
        eSTAT = 'x';
    end 
    
    %-Retrieve design matrix.
    %----------------------------------------------------------------------
    X = SwE.xX.X;
    
    %-Create the xCon object.
    %----------------------------------------------------------------------
    DxCon = spm_FcUtil('Set', 'Contrast 1: Activation', STAT, 'c', c', X);
    
    % Work out if we are in clusterwise bootstrap or not.
    %----------------------------------------------------------------------
    if SwE.WB.clusterWise
        wbstring = '-WB';
    else
        wbstring = '';
    end
    
    % Add the SwE volumes.
    %----------------------------------------------------------------------
    DxCon.Vspm = spm_vol(['swe_vox_' eSTAT STAT 'stat_c01' file_ext]);
%    DxCon.Vspm2 = spm_vol(['swe_vox_' STAT 'stat' wbstring '_c01' file_ext]); %THIS DOESN'T EXIST FOR WB - MUST FIX ASAP
    DxCon.VspmUncP = spm_vol(['swe_vox_' STAT 'stat_lp' wbstring '_c01' file_ext]);
    DxCon.VspmFDRP = spm_vol(['swe_vox_' STAT 'stat_lpFDR' wbstring '_c01' file_ext]);
    DxCon.VspmFWEP = spm_vol(['swe_vox_' STAT 'stat_lpFWE' wbstring '_c01' file_ext]);
    DxCon.VspmFWEP_clus = spm_vol(['swe_clustere_' STAT 'stat_lpFWE' wbstring '_c01' file_ext]);
    
    % Return SwE.
    %----------------------------------------------------------------------
    if ~isfield(SwE, 'xCon')
        SwE.xCon = DxCon;
    else
        SwE.xCon = [SwE.xCon DxCon];
    end

end