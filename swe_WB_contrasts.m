% This function creates the missing 'xCon' needed for wild bootstrap SwE
% objects.
%
% Author: Tom Maullin (08/06/2018)
function [SwE] = swe_WB_contrasts(SwE)

    %-Get and change to results directory
    %----------------------------------------------------------------------
    try
        cd(SwE.swd);
    end
    
    % Retrieve contrast details.
    STAT = SwE.WB.stat;
    c = SwE.WB.con;
    
    % Retrieve design matrix.
    X = SwE.xX.X;
    
    % Create the xCon object.
    DxCon = spm_FcUtil('Set', 'Contrast 1: Activation', STAT, 'c', c', X);
    disp(DxCon);
    
    % Add the SwE volumes.
    DxCon.Vedf = 
    

end