function [SwE] = swe_contrasts_WB(SwE)
% This function creates the missing 'xCon' for wild bootstrap SwE objects.
% =========================================================================
% Note: Whilst contrasts are created in 'swe_contrasts.m', in
% 'swe_WB_contrasts.m' they are only recorded. They have already been
% created by 'swe_cp_WB.m'.
% =========================================================================
% FORMAT [SwE] = SwE_contrasts_WB(SwE)
% -------------------------------------------------------------------------
% Inputs/Outputs:
%   - SwE: SwE data structure
% =========================================================================
% Author: Tom Maullin (08/06/2018)
% Version Info:  $Format:%ci$ $Format:%h$

    %-Get and change to results directory
    %----------------------------------------------------------------------
    try
        cd(SwE.swd);
    end

    %-Get file extension
    %----------------------------------------------------------------------
    file_ext = swe_get_file_extension(SwE.xY.P{1});
    isMat          = strcmpi(file_ext,'.mat');
    isCifti  = strcmpi(file_ext,'.dtseries.nii') ||  strcmpi(file_ext,'.dscalar.nii');

    if isCifti
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
        else
            file_ext = spm_file_ext;
            file_data_type = 'vox';
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
    if isfield(SwE,'WB')
        wbstring = '-WB';
    else
        wbstring = '';
    end

    % Add the SwE volumes.
    %----------------------------------------------------------------------
    DxCon.Vspm = swe_data_hdr_read(sprintf('swe_%s_%c%cstat_c%.2d%s', file_data_type, eSTAT, STAT, 1, file_ext));
    DxCon.VspmUncP = swe_data_hdr_read(sprintf('swe_%s_%cstat_lp%s_c%.2d%s', file_data_type, STAT, wbstring, 1, file_ext));
    DxCon.VspmFDRP = swe_data_hdr_read(sprintf('swe_%s_%cstat_lpFDR%s_c%.2d%s', file_data_type, STAT, wbstring, 1, file_ext));
    DxCon.VspmFWEP = swe_data_hdr_read(sprintf('swe_%s_%cstat_lpFWE%s_c%.2d%s', file_data_type, STAT, wbstring, 1, file_ext));
    if isfield(SwE,'WB')
        if SwE.WB.clusterWise
            % In the most recent versions of the toolbox, VspmFWEP_clusnorm should be used, but, for models estimated
            % with older versions, it does not exist and VspmFWEP_clus should be used instead.
            try
                DxCon.VspmFWEP_clusnorm = swe_data_hdr_read(sprintf('swe_clusternorm_%cstat_lpFWE%s_c%.2d%s', STAT, wbstring, 1, file_ext));
            catch
                DxCon.VspmFWEP_clus = swe_data_hdr_read(sprintf('swe_clustere_%cstat_lpFWE%s_c%.2d%s', STAT, wbstring, 1, file_ext));
            end
        end
        if isfield(SwE.WB, 'TFCE')
            DxCon.VspmTFCE = swe_data_hdr_read(sprintf('swe_tfce_c%.2d%s', 1, file_ext));
            DxCon.VspmTFCEP = swe_data_hdr_read(sprintf('swe_tfce_lp-WB_c%.2d%s', 1, file_ext));
            DxCon.VspmTFCEFWEP = swe_data_hdr_read(sprintf('swe_tfce_lpFWE-WB_c%.2d%s', 1, file_ext));
        end
    end
    DxCon.Vedf = swe_data_hdr_read(sprintf('swe_%s_edf_c%.2d%s', file_data_type, 1, file_ext));

    % Add the positive contrast.
    %----------------------------------------------------------------------
    if ~isfield(SwE, 'xCon')
        SwE.xCon = DxCon;
    else
        SwE.xCon = [SwE.xCon DxCon];
    end

    % We need to add the negative contrast if we are doing a T test.
    if strcmpi(STAT, 't')

        DxCon = spm_FcUtil('Set', 'Contrast 2: Deactivation', STAT, 'c', -c', X);

        DxCon.Vspm = swe_data_hdr_read(sprintf('swe_%s_%c%cstat_c%.2d%s', file_data_type, eSTAT, STAT, 2, file_ext));
        DxCon.VspmUncP = swe_data_hdr_read(sprintf('swe_%s_%cstat_lp%s_c%.2d%s', file_data_type, STAT, wbstring, 2, file_ext));
        DxCon.VspmFDRP = swe_data_hdr_read(sprintf('swe_%s_%cstat_lpFDR%s_c%.2d%s', file_data_type, STAT, wbstring, 2, file_ext));
        DxCon.VspmFWEP = swe_data_hdr_read(sprintf('swe_%s_%cstat_lpFWE%s_c%.2d%s', file_data_type, STAT, wbstring, 2, file_ext));
        if isfield(SwE,'WB')
            if SwE.WB.clusterWise
                % In the most recent versions of the toolbox, VspmFWEP_clusnorm should be used, but, for models estimated
                % with older versions, it does not exist and VspmFWEP_clus should be used instead.
                try
                    DxCon.VspmFWEP_clusnorm = swe_data_hdr_read(sprintf('swe_clusternorm_%cstat_lpFWE%s_c%.2d%s', STAT, wbstring, 2, file_ext));
                catch
                    DxCon.VspmFWEP_clus = swe_data_hdr_read(sprintf('swe_clustere_%cstat_lpFWE%s_c%.2d%s', STAT, wbstring, 2, file_ext));
                end
            end
            if isfield(SwE.WB, 'TFCE')
                DxCon.VspmTFCE = swe_data_hdr_read(sprintf('swe_tfce_c%.2d%s', 2, file_ext));
                DxCon.VspmTFCEP = swe_data_hdr_read(sprintf('swe_tfce_lp-WB_c%.2d%s', 2, file_ext));
                DxCon.VspmTFCEFWEP = swe_data_hdr_read(sprintf('swe_tfce_lpFWE-WB_c%.2d%s', 2, file_ext));
            end
        end
        DxCon.Vedf = swe_data_hdr_read(sprintf('swe_%s_edf_c%.2d%s', file_data_type, 1, file_ext));

        SwE.xCon = [SwE.xCon DxCon];
    end

end
