function varargout = swe(varargin)
% Sandwich Estimator for Neuroimaging Longitudinal Data Analysis, SwE.
%
% This function initializes things for swe and provides some low level
% functionalities

% Written by Bryan Guillaume
if nargin == 0,
    Action = 'StartUp';
else
    Action = varargin{1};
end

switch lower(Action)
    %==================================================================
    case 'startup'                                % Startup the toolbox
    %==================================================================     
        % check installation of toolbox and of SPM8/SPM12
        ok = check_installation;
        if ~ok
            beep
            fprintf('INSTALLATION PROBLEM!');
            return
        end
        
        % Welcome message
        swe('ASCIIwelcome');
        
        % launch the main GUI
        swe_ui_main;
        
        % print present working directory
        fprintf('SwE present working directory:\n\t%s\n',pwd)
        
        
        %==================================================================
    case 'asciiwelcome'                          %-ASCII swe banner welcome
        %==================================================================
        disp( '   ___          ___           ___     _                   ');
        disp( '  / __) _    _ | __)   _  _  | _ |   /_|                  ');
        disp( '  \__ \ \\/\// | __)   \\//  ||_|| _  ||                  ');
        disp( '  (___/  \/\/  |___)    \/   |___||_| ||                  ');
        fprintf('\n  swe v0.1\n');
  
        
        %==================================================================
    otherwise                                       %-Unknown action string
        %==================================================================
        error('Unknown action string');
        
end

return

%=======================================================================
%% SUBFUNCTIONS
%=======================================================================

%=======================================================================
function ok = check_installation
% function to check installation state of toolbox,
% particullarly the SPM path setup

ok=1;

if exist('spm.m','file')
    if ~strcmpi(spm('ver'),'spm8')&& ~strncmpi(spm('ver'),'spm12',5)
        beep
        fprintf('\nERROR:\n')
        fprintf('\tSPM8 or SPM12 should be installed on your computer, and\n')
        fprintf('\tbe available on MATLABPATH!\n\n')
        ok = 0;
    end
else
    beep
    fprintf('\nERROR:\n')
    fprintf('\tSPM8 or SPM12 should be installed on your computer, and\n')
    fprintf ('\tbe available on MATLABPATH!\n\n')
    ok = 0;
end

return