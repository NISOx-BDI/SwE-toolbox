function varargout = swe_text_input(varargin)
% swe_TEXT_INPUT M-file for swe_text_input.fig
% 
% swe_TEXT_INPUT, by itself, creates a new swe_TEXT_INPUT or raises the 
% existing singleton*.
%
% H = swe_TEXT_INPUT returns the handle to a new swe_TEXT_INPUT or the 
% handle to the existing singleton*.
%
% swe_TEXT_INPUT('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in swe_TEXT_INPUT.M with the given input arguments.
%
% swe_TEXT_INPUT('Property','Value',...) creates a new swe_TEXT_INPUT or 
% raises the existing singleton*.  Starting from the left, property value 
% pairs are applied to the GUI before swe_text_input_OpeningFcn gets called.
% An unrecognized property name or invalid value makes property application
% stop.  All inputs are passed to swe_text_input_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Written by Bryan Guillaume

% Edit the above text to modify the response to help swe_text_input

% Last Modified by GUIDE v2.5 06-Apr-2013 03:53:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @swe_text_input_OpeningFcn, ...
                   'gui_OutputFcn',  @swe_text_input_OutputFcn, ...
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

% --- Executes just before swe_text_input is made visible.
function swe_text_input_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to swe_text_input (see VARARGIN)

%if window already exists, just put it as the current figure
Tag='text_input';
F = findall(allchild(0),'Flat','Tag',Tag);
if length(F) > 1
    % Multiple Graphics windows - close all but most recent
    close(F(2:end))
    F = F(1);
    uistack(F,'top')
elseif length(F)==1
    uistack(F,'top')
else
    set(handles.figure1,'Tag',Tag)
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
    x=get(handles.figure1,'Position');
    set(handles.figure1,'DefaultTextFontSize',FS*12,...
        'DefaultUicontrolFontSize',FS*12,...
        'DefaultTextFontName',PF,...
        'DefaultAxesFontName',PF,...
        'DefaultUicontrolFontName',PF)
    set(handles.figure1,'Position',ratio.*x)
    set(handles.figure1,'Resize','on')
    
    % Choose some figure parameters
    aa=get(handles.figure1,'children');
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
        set(aa(i),'FontUnits','pixel')
        xf=get(aa(i),'FontSize');
        if ispc
            set(aa(i),'FontSize',ceil(FS*xf),'FontName',PF,...
                'FontUnits','normalized','Units','normalized')
        else
            set(aa(i),'FontSize',ceil(FS*xf),'FontName',PF,...
                'Units','normalized')
        end
    end
    
    
    % Choose default command line output for swe_text_input
    handles.output = hObject;
    % Update handles structure
    guidata(hObject, handles);
    if ~isempty(varargin) && strcmpi(varargin{1},'Title')
        set(gcf,'Name',varargin{2});
        if length(varargin)>3 && strcmpi(varargin{3},'UserData')
            set(handles.edit1,'String',varargin{4})
            handles.output = get(handles.edit1,'String');
        end
    end
end
% Update handles structure
guidata(hObject, handles);

%UIWAIT makes swe_text_input wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = swe_text_input_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isfield(handles,'output') && ~isempty(handles.output)
    varargout{1} = handles.output;
else
    varargout{1}=[];
end

% The figure can be deleted now
if isfield(handles,'figure1')
    delete(handles.figure1);
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in okbutton.
function okbutton_Callback(hObject, eventdata, handles)
% hObject    handle to okbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);

% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
handles.output=[];
% Update handles structure
guidata(hObject, handles);

uiresume(handles.figure1);
