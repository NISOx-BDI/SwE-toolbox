% This function creates the missing 'xCon' needed for wild bootstrap SwE
% objects. Note: Whilst contrasts are created in 'swe_contrasts.m', in 
% 'swe_WB_contrasts.m' they are only recorded. They have already been
% created by 'swe_cp_WB.m'.
% =========================================================================
% FORMAT [SwE] = SwE_contrasts_WB(SwE)
% -------------------------------------------------------------------------
% Inputs/Outputs:
%   - SwE: SwE data structure
% =========================================================================
% Author: Tom Maullin (08/06/2018)

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
    DxCon.Vspm = spm_vol(sprintf('swe_vox_%c%cstat_c%.2d%s', eSTAT, STAT, 1, file_ext));
    DxCon.VspmUncP = spm_vol(sprintf('swe_vox_%cstat_lp_c%.2d%s', STAT, 1, file_ext));
    DxCon.VspmFDRP = spm_vol(sprintf('swe_vox_%cstat_lpFDR%s_c%.2d%s', STAT, wbstring, 1, file_ext));
    DxCon.VspmFWEP = spm_vol(sprintf('swe_vox_%cstat_lpFWE%s_c%.2d%s', STAT, wbstring, 1, file_ext));
    DxCon.VspmFWEP_clus = spm_vol(sprintf('swe_clustere_%cstat_lpFWE%s_c%.2d%s', STAT, wbstring, 1, file_ext));
    DxCon.Vedf = spm_vol(sprintf('swe_vox_edf_c%.2d%s', 1, file_ext));
    
    % Return SwE.
    %----------------------------------------------------------------------
    if ~isfield(SwE, 'xCon')
        SwE.xCon = DxCon;
    else
        SwE.xCon = [SwE.xCon DxCon]; 
    end

end