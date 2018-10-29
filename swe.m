function varargout = swe(varargin)
% This function initializes the SwE toolbox and checks the installation.
% =========================================================================
% FORMAT swe, swe(), swe('startup') - starts the SwE toolbox.
% FORMAT swe('ASCIIwelcome') - display the startup message.
% FORMAT swe('Colour') - returns the interface colour.
% -------------------------------------------------------------------------
% Inputs:
%  - str - 'ASCIIwelcome' for startup message or 'startup' for toolbox
% =========================================================================
% Written by Bryan Guillaume
% $Id$

versionNo = '2.0.0.rc';

if nargin == 0
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
        
        % Add pathes for SPM functions
        addpath(fullfile(spm('Dir'),'matlabbatch'))
        if strncmpi(spm('ver'),'spm12',5)
            addpath(fullfile(spm('Dir'),'compat'))
        end
        
        % Add path to SwE toolbox.
        addpath(fileparts(mfilename('fullpath')));
        
        % launch the main GUI
        swe_ui_main;
        
        % print present working directory
        fprintf('SwE present working directory:\n\t%s\n',pwd)
        
        
        %==================================================================
    case 'asciiwelcome'                          %-ASCII swe banner welcome
        %==================================================================
        disp( '   ___          ___           _      ___     _  _        ');
        disp( '  / __) _    _ | __)   _  _  /_|    (__ \   /_|/_|       ');
        disp( '  \__ \ \\/\// | __)   \\//   || _  / __/ _  || ||       ');
        disp( '  (___/  \/\/  |___)    \/    |||_| \___)|_| || ||       ');
        fprintf('\n  swe v%s \n', versionNo);
  
    case 'colour'
        
        varargout = {[0.8 0.8 1.0], 'Diluted Blackcurrent Purple'};
        
    case 'ver'
        
        varargout{1}=versionNo;
        
        %==================================================================
    otherwise                                       %-Unknown action string
        %==================================================================
        error('Unknown action string');
        
end

return

%=======================================================================
% SUBFUNCTIONS
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