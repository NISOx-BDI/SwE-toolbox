function varargout = swe_ui_main(varargin)
% SWE_UI_MAIN MATLAB code for SWE_UI_MAIN.fig
%      SWE_UI_MAIN, by itself, creates a new SWE_UI_MAIN or raises the existing
%      singleton*.
%
%      H = SWE_UI_MAIN returns the handle to a new SWE_UI_MAIN or the handle to
%      the existing singleton*.
%
%      SWE_UI_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SWE_UI_MAIN.M with the given input arguments.
%
%      SWE_UI_MAIN('Property','Value',...) creates a new SWE_UI_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SWE_UI_MAIN_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SWE_UI_MAIN_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Written by Bryan Guillaume

% Last Modified by GUIDE v2.5 24-Apr-2013 18:41:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @swe_ui_main_OpeningFcn, ...
                   'gui_OutputFcn',  @swe_ui_main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before swe_ui_main is made visible.
function swe_ui_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to swe_ui_main (see VARARGIN)

Tag='swemain';
F = findall(allchild(0),'Flat','Tag',Tag);
if length(F) > 1
    % Multiple Graphics windows - close all but most recent
    close(F(2:end))
    F = F(1);
    uistack(F,'top')
elseif length(F)==1
    uistack(F,'top')
else
    
    %set size of the window, taking screen resolution and platform into account
    S0= spm('WinSize','0',1);   %-Screen size (of the current monitor)
    if ispc
        PF='MS Sans Serif';
    else
        PF= spm_platform('fonts');     %-Font names (for this platform)
        PF=PF.helvetica;
    end
    tmp  = [S0(3)/1280 (S0(4))/800];
    ratio=min(tmp)*[1 1 1 1];
    FS = 1 + 0.85*(min(ratio)-1);  %factor to scale the fonts
    x=get(handles.swemain,'Position');
    set(handles.swemain,'DefaultTextFontSize',FS*8,...
        'DefaultUicontrolFontSize',FS*8,...
        'DefaultTextFontName',PF,...
        'DefaultAxesFontName',PF,...
        'DefaultUicontrolFontName',PF)
    set(handles.swemain,'Position',ratio.*x)
    
    % Choose some figure parameters
    aa=get(handles.swemain,'children');
    for i=1:length(aa)
        if strcmpi(get(aa(i),'type'),'uipanel')
            bb=get(aa(i),'children');
            if ~isempty(bb)
                for j=1:length(bb)
                    set(bb(j),'FontUnits','pixel')
                    xf=get(bb(j),'FontSize');
                    set(bb(j),'FontSize',ceil(FS*xf),'FontName',PF,...
                        'FontUnits','normalized','Units','normalized')
                end
            end
        end
        xf=get(aa(i),'FontSize');
        if ispc
            set(aa(i),'FontSize',FS*xf,'FontName',PF,...
                'FontUnits','normalized','Units','normalized')
        else
            set(aa(i),'FontSize',FS*xf,'FontName',PF,...
                'Units','normalized')
        end
    end    
end


% Choose default command line output for swe_ui_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes swe_ui_main wait for user response (see UIRESUME)
% uiwait(handles.swemain);


% --- Outputs from this function are returned to the command line.
function varargout = swe_ui_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in specify.
function specify_Callback(hObject, eventdata, handles)
% hObject    handle to specify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
swe_design



% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
swe_cp


% --- Executes on button press in results.
function results_Callback(hObject, eventdata, handles)
% hObject    handle to results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
swe_results_ui;


% --- Executes on button press in batch.
function batch_Callback(hObject, eventdata, handles)
% hObject    handle to batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
swe_batch
