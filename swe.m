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
% Version Info:  $Format:%ci$ $Format:%h$

versionNo = '2.0.0';

try
  Modality = spm_get_defaults('modality');
catch
  spm('defaults','FMRI'); % Starting w/out SPM
end

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
        a = generateAscii(['SwE v' versionNo]);
        fprintf('%s \n', a{1}, a{2}, a{3}, a{4});
        fprintf('swe v%s \n', versionNo);
  
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

end

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

end

% The following functions are for generating the ascii welcome message for
% SwE. 
% -------------------------------------------------------------------------
% This method converts character 'char' to it's equivalent 4-line ascii
% art. New characters can be added by creating new cases in the below
% switch.
function aChar = char2ascii(char)
    
    switch char
        
        case '0'
            
            aChar = {' ___ ',...
                     '|   |',...
                     '| | |',...
                     '|___|'};
            
        case '1'
            
            aChar = {' _ ',...
                     '/_|',...
                     ' ||',...
                     ' ||'};
                
        case '2'
            
            aChar = {' ___ ',...
                     '(__ \',...
                     '/ __/',...
                     '\___)'};
                
        case '3'
            
            aChar = {' ___ ',...
                     '(__ )',...
                     ' (_ \',...
                     '(___/'};           
        case '4'
            
            aChar = {'  __ ',...
                     ' /. |',...
                     '(_  _)',...
                     '  (_) '};           
        case '5' 
            
            aChar = {' ___ ',...
                     '| __)',...
                     '|__ \',...
                     '(___/'};
                
        case '6'
            
            aChar = {'  _  ',...
                     ' / ) ',...
                     '/ _ \',...
                     '\___/'};     
        case '7'
            
            aChar = {' ___ ',...
                     '(__ )',...
                     ' / / ',...
                     '(_/  '};      
        case '8'
            
            aChar = {' ___ ',...
                     '| _ |',...
                     '| _ |',...
                     '|___| '};   
        case '9'
            
            aChar = {' ___ ',...
                     '/ _ \',...
                     '\_  /',...
                     ' (_/ '};   
                
        case '.'
            
            aChar = {'   ',...
                     '   ',...
                     ' _ ',...
                     '|_|'};  
                
        case 'S'
 
            aChar = {' ___ ',...
                     '/ __)',...
                     '\__ \',...
                     '(___/'}; 
                
         case 'w'
 
            aChar = {'      ',...
                     '_    _',...
                     '\\/\//',...
                     ' \/\/ '};  
                
        case 'E'
 
            aChar = {' ___ ',...
                     '| __)',...
                     '| __)',...
                     '|___)'};      
                
        case ' '
            
            aChar = {'   ',...
                     '   ',...
                     '   ',...
                     '   '};
                
        case 'v'
            
            aChar = {'     ',...
                     '_  _ ',...
                     '\\// ',...
                     ' \/  '};  
            
            
    end
    
end

% This function takes two cell arrays containing ascii art in the form
% {line1, line2, line3,...} and concatenates them horizontally.
function c = concatAscii(a,b)

    c = arrayfun(@(i) [a{i} b{i}], 1:length(b),'UniformOutput',false);
end

% This function takes in a string as input and recursively generates the
% ASCII art representation of said string.
function ascii = generateAscii(str)
    
    if length(str)~=1
        ascii = concatAscii(generateAscii(str(1:(end-1))), generateAscii(str(end)));
    else
        ascii = char2ascii(str);
    end
    
end