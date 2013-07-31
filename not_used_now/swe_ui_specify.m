function varargout = swe_ui_specify(varargin)
% swe_UI_SPECIFY M-file for swe_ui_specify.fig
% 
% swe_UI_SPECIFY, by itself, creates a new swe_UI_SPECIFY or raises the 
% existing singleton*.
%
% H = swe_UI_SPECIFY returns the handle to a new swe_UI_SPECIFY or the handle
% to the existing singleton*.
%
% swe_UI_SPECIFY('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in swe_UI_SPECIFY.M with the given input arguments.
%
% swe_UI_SPECIFY('Property','Value',...) creates a new swe_UI_SPECIFY or 
% raises the existing singleton*.  Starting from the left, property value 
% pairs are applied to the GUI before swe_ui_specify_OpeningFcn gets called.
% An unrecognized property name or invalid value makes property application
% stop.  All inputs are passed to swe_ui_specify_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________


% Edit the above text to modify the response to help swe_ui_specify

% Last Modified by GUIDE v2.5 10-Apr-2013 16:08:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @swe_ui_specify_OpeningFcn, ...
                   'gui_OutputFcn',  @swe_ui_specify_OutputFcn, ...
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


% --- Executes just before swe_ui_specify is made visible.
function swe_ui_specify_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to swe_ui_specify (see VARARGIN)

% Choose default command line output for swe_ui_specify
handles.output = hObject;
%if window already exists, just put it as the current figure
Tag='swe_specify';
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
        xf=get(aa(i),'FontSize');
        if ispc
            set(aa(i),'FontSize',ceil(FS*xf),'FontName',PF,...
                'FontUnits','normalized','Units','normalized')
        else
            set(aa(i),'FontSize',ceil(FS*xf),'FontName',PF,...
                'Units','normalized')
        end
    end
    handles.saved=0;
    handles.cgr=1; %current group
    handles.cs=1;
    handles.cv=1;
    handles.cf=1;
    handles.dat=struct('dir',[],'group',[],'design',[],'mask',[]);
    handles.dat.mask.val=0; % set flag "no mask"
    handles.dat.mask.fname={}; % set flag "no mask"
    handles.dat.swe_type='Hom.';
    handles.dat.SS=0;  

end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes swe_ui_specify wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = swe_ui_specify_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in group_list.
function group_list_Callback(hObject, eventdata, handles)
% hObject    handle to group_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns group_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from group_list
handles.cgr=get(handles.group_list,'Value');
handles.cs=1;
handles.cv=1;
handles.cf=1;
update_display_data(hObject,handles);

handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function group_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to group_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in subjects_list.
function subjects_list_Callback(hObject, eventdata, handles)
% hObject    handle to subjects_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns subjects_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        subjects_list
handles.cs=get(handles.subjects_list,'Value');
handles.cv=1;
handles.cf=1;
update_display_data(hObject,handles);
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function subjects_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subjects_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in visits_list.
function visits_list_Callback(hObject, eventdata, handles)
% hObject    handle to visits_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns visits_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from visits_list
handles.cv=get(handles.visits_list,'Value');
handles.cf=1;
update_display_data(hObject,handles);
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function visits_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visits_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in files_list.
function files_list_Callback(hObject, eventdata, handles)
% hObject    handle to files_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns files_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        files_list
handles.cf=get(handles.files_list,'Value');
update_display_data(hObject,handles);
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function files_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to files_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in br_res_dir.
function br_res_dir_Callback(hObject, eventdata, handles)
% hObject    handle to br_res_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dat.dir=uigetdir(cd,'Directory to write results');
set(handles.edit1,'String',handles.dat.dir,'FontAngle','normal')
handles.saved=0;
% Update handles structure
guidata(hObject, handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.dat.dir=get(handles.edit1,'String');
handles.saved=0;
% Update handles structure
guidata(hObject, handles);

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

% --- Executes on button press in gr_add.
function gr_add_Callback(hObject, eventdata, handles)
% hObject    handle to gr_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'ds')
    defnam=['G',num2str(length(handles.ds)+1)];
    it=1;
    namlist={handles.dat.group(:).gr_name};
    while any(strcmpi(namlist,defnam))
       defnam=['G',num2str(length(handles.ds)+1+it)]; 
       it=it+1;
    end
else
    defnam='G1';
end
gname=swe_text_input('Title','Enter group name','UserData',defnam);
if isnumeric(gname)
    return
end
if isfield(handles.dat,'group') && ...
        isfield(handles.dat.group,'gr_name')
    namlist={handles.dat.group(:).gr_name};
    if any(strcmpi(namlist,gname))
        it=1;
        while any(strcmpi(namlist,gname))
            gname=['G', num2str(length(namlist)+it)];
            it=it+1;
        end 
        disp(['This group has already been defined, name attributed is ', ...
            'G', num2str(length(namlist)+1)]);
        beep
    end
end
if ~isfield(handles,'ds')
    handles.ds=cell(1);
else
    handles.ds=[handles.ds; cell(1)];
end
ngr=length(handles.ds);
handles.cgr=ngr;
handles.dat.group(ngr).gr_name=gname;
newlist=[get(handles.group_list,'String'); {gname}];
set(handles.group_list,'String',newlist);
ren=uicontextmenu;
% Update handles structure
guidata(hObject, handles);
item1=uimenu(ren,'Label','Rename','Callback',@rengroup);
set(handles.group_list,'UIContextMenu',ren)
handles=guidata(hObject);
set(handles.subj_add,'enable','on')
set(handles.subj_remove,'enable','on')

% Update handles structure
guidata(hObject, handles);

update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
% Update handles structure
guidata(hObject, handles);


%Function called when right-clicking on the 'rename' menu
function rengroup(hObject,eventdata)
handles=guidata(hObject);
val=get(handles.group_list,'Value');
renam=swe_text_input('Title','Rename group');
if isempty(renam)
    return
end
handles.dat.group(val).gr_name=renam;
list=get(handles.group_list,'String');
list{val}=renam;
set(handles.group_list,'String',list)
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in gr_remove.
function gr_remove_Callback(hObject, eventdata, handles)
% hObject    handle to gr_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list=get(handles.group_list,'String');
ngr=length(handles.ds);
if ngr==1
    nlist=[];
    handles.dat.group=[];
    handles.cs=0;
    handles.cv=0;
    handles.cf=0;
    handles=rmfield(handles,'ds');
    handles.visits_list={};
elseif handles.cgr==1 && ngr>1
    nlist=list(2:end);
    handles.dat.group=handles.dat.group(2:end);
    handles.ds=handles.ds(2:end);
    handles.cgr=1;
    handles.cs=1;
    handles.cv=1;
    handles.cf=1;
elseif handles.cgr==ngr
    nlist=list(1:end-1);
    handles.dat.group=handles.dat.group(1:end-1);
    handles.cgr=handles.cgr-1;
    handles.ds=handles.ds(1:end-1);
    handles.cs=1;
    handles.cv=1;
    handles.cf=1;
else
    nlist=list([1:handles.cgr-1,handles.cgr+1:end]);
    handles.ds=handles.ds([1:handles.cgr-1,handles.cgr+1:end]);
    handles.dat.group=handles.dat.group([1:handles.cgr-1,handles.cgr+1:end]);
    handles.cgr=handles.cgr-1;
    handles.cs=1;
    handles.cv=1;
    handles.cf=1;
end
set(handles.group_list,'String',nlist);

set(handles.subj_add,'enable','on')
set(handles.subj_remove,'enable','on')

update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in subj_add.
function subj_add_Callback(hObject, eventdata, handles)
% hObject    handle to subj_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'ds') || length(handles.ds)<1
    beep
    disp('Please add at least one group before adding subjects')
    return
end
defnam=['S',num2str(length(handles.ds{handles.cgr})+1)];
if ~isempty(handles.ds{handles.cgr})
    it=1;
    namlist={handles.dat.group(handles.cgr).subject(:).subj_name};
    while any(strcmpi(namlist,defnam))
        defnam=['S',num2str(length(handles.ds{handles.cgr})+1+it)];
        it=it+1;
    end
end
sname=swe_text_input('Title','Enter subject name','UserData',defnam);
if isnumeric(sname)
    return
elseif ~isempty(strfind(lower(sname),'scan'))
    beep
    disp('Scan(s) is a reserved name. Please correct')
    return
end
if isfield(handles.dat.group(handles.cgr),'subject') && ...
        isfield(handles.dat.group(handles.cgr).subject,'subj_name')
    namlist={handles.dat.group(handles.cgr).subject(:).subj_name};
    if any(strcmpi(namlist,sname))
        it=1;
        while any(strcmpi(namlist,sname))
            sname=['S', num2str(length(namlist)+it)];
            it=it+1;
        end 
        disp(['This subject has already been defined, name attributed is ', ...
            'S', num2str(length(namlist)+1)]); 
        beep
    end
end
handles.ds{handles.cgr}=[handles.ds{handles.cgr}, cell(1)];
handles.cs=length(handles.ds{handles.cgr});
if ~isfield(handles.dat.group(handles.cgr),'subject')
    handles.dat.group(handles.cgr).subject=[];
end
handles.dat.group(handles.cgr).subject(handles.cs).subj_name=sname;
newlist=[get(handles.subjects_list,'String'); {sname}];
set(handles.subjects_list,'String',newlist);
rens=uicontextmenu;
% Update handles structure
guidata(hObject, handles);
item1=uimenu(rens,'Label','Rename','Callback',@rensubj);
set(handles.subjects_list,'UIContextMenu',rens)
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);
update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
% Update handles structure
guidata(hObject, handles);

%Function called when right-clicking on the 'rename' menu
function rensubj(hObject,eventdata)
handles=guidata(hObject);
val=get(handles.subjects_list,'Value');
renam=swe_text_input('Title','Rename subject');
if isempty(renam)
    return
end
if ~isempty(strfind(lower(renam),'scan'))
    beep
    disp('Scan(s) is a reserved name. Please correct')
    return
end
handles.dat.group(handles.cgr).subject(val).subj_name=renam;
list=get(handles.subjects_list,'String');
list{val}=renam;
set(handles.subjects_list,'String',list)
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in subj_remove.
function subj_remove_Callback(hObject, eventdata, handles)
% hObject    handle to subj_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list=get(handles.subjects_list,'String');
cgr=handles.cgr;
nsubj=length(handles.ds{cgr});
if nsubj==1
    nlist=[];
    handles.dat.group(cgr).subject=[];
    handles.cs=0;
    handles.cv=0;
    handles.cf=0;
    handles.ds(cgr)={};
    handles.visits_list={};
elseif handles.cs==1 && nsubj>1
    nlist=list(2:end);
    handles.dat.group(cgr).subject=handles.dat.group(cgr).subject(2:end);
    handles.ds{cgr}=handles.ds{cgr}(2:end);
    handles.cs=1;
    handles.cv=1;
    handles.cf=1;
elseif handles.cs==nsubj
    nlist=list(1:end-1);
    handles.dat.group(cgr).subject=handles.dat.group(cgr).subject(1:end-1);
    handles.ds{cgr}=handles.ds{cgr}(1:end-1);
    handles.cs=handles.cs-1;
    handles.cv=1;
    handles.cf=1;
else
    nlist=list([1:handles.cs-1,handles.cs+1:end]);
    handles.dat.group(cgr).subject=handles.dat.group(cgr).subject([1:handles.cs-1,handles.cs+1:end]);
    handles.ds{cgr}=handles.ds{cgr}([1:handles.cs-1,handles.cs+1:end]);
    handles.cs=handles.cs-1;
    handles.cv=1;
    handles.cf=1;
end
set(handles.subjects_list,'String',nlist);
update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in vis_add.
function vis_add_Callback(hObject, eventdata, handles)
% hObject    handle to vis_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'ds') || length(handles.ds)<1 || ...
        length(handles.ds(handles.cgr))<1
    beep
    disp(['Please add at least one group and one subject before adding visits'])
    return
end
defnam=['V',num2str(length(handles.ds{handles.cgr}{handles.cs})+1)];
if ~isempty(handles.ds{handles.cgr}{handles.cs})
    it=1;
    namlist={handles.dat.group(handles.cgr).subject(handles.cs).visit(:).vis_name};
    while any(strcmpi(namlist,defnam))
        defnam=['V',num2str(length(handles.ds{handles.cgr}{handles.cs})+1+it)];
        it=it+1;
    end
end
vname=swe_text_input('Title','Enter visit name','UserData',defnam);
if isnumeric(vname)
    return
elseif ~isempty(strfind(lower(vname),'scan'))
    beep
    disp('Scan(s) is a reserved name. Please correct')
    return
end
if isfield(handles.dat.group(handles.cgr),'subject') && ...
        isfield(handles.dat.group(handles.cgr).subject(handles.cs),'visit') && ...
            isfield(handles.dat.group(handles.cgr).subject(handles.cs).visit,'vis_name')
    namlist={handles.dat.group(handles.cgr).subject(handles.cs).visit(:).vis_name};
    if any(strcmpi(namlist,vname))
        it=1;
        while any(strcmpi(namlist,vname))
            vname=['V', num2str(length(namlist)+it)];
            it=it+1;
        end 
        disp(['This visit has already been defined, name attributed is ', ...
            'V', num2str(length(namlist)+1)]);
        beep
    end
end
handles.ds{handles.cgr}{handles.cs}=[handles.ds{handles.cgr}{handles.cs}; cell(1)];
handles.cv=length(handles.ds{handles.cgr}{handles.cs});
if ~isfield(handles.dat.group(handles.cgr).subject(handles.cs),'visit')
    handles.dat.group(handles.cgr).subject(handles.cs).visit=[];
end
handles.dat.group(handles.cgr).subject(handles.cs).visit(handles.cv).vis_name=vname;
newlist=[get(handles.visits_list,'String'); {vname}];
set(handles.visits_list,'String',newlist);
renv=uicontextmenu;
% Update handles structure
guidata(hObject, handles);
item1=uimenu(renv,'Label','Rename','Callback',@renvis);
set(handles.visits_list,'UIContextMenu',renv)
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);

update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
% Update handles structure
guidata(hObject, handles);


%Function called when right-clicking on the 'rename' menu
function renvis(hObject,eventdata)
handles=guidata(hObject);
val=get(handles.visits_list,'Value');
renam=swe_text_input('Title','Rename visit');
if isempty(renam)
    return
end
if ~isempty(strfind(lower(renam),'scan'))
    beep
    disp('Scan(s) is a reserved name. Please correct')
    return
end
handles.dat.group(handles.cgr).subject(handles.cs).visit(val).vis_name=renam;
list=get(handles.visits_list,'String');
list{val}=renam;
set(handles.visits_list,'String',list)
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in vis_remove.
function vis_remove_Callback(hObject, eventdata, handles)
% hObject    handle to vis_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vlist=get(handles.visits_list,'String');
if ~isempty(get(handles.files_list,'String'))
    flist=get(handles.files_list,'String');
end
cgr=handles.cgr;
cs=handles.cs;
cv=handles.cv;
nvis=length(handles.ds{cgr}{cs});
if nvis==1
    vlist=[];
    handles.dat.group(cgr).subject(cs).visit=[];
    handles.cv=0;
    handles.cf=0;
    handles.ds{cgr}{cs}={};
    if ~isempty(get(handles.files_list,'String'))
        flist={};
        handles.dat.group(cgr).subject(cs).files=[];
    end
elseif cv==1 && nvis>1
    vlist=vlist(2:end);
    handles.dat.group(cgr).subject(cs).visit=handles.dat.group(cgr).subject(cs).visit(2:end);
    handles.ds{cgr}{cs}=handles.ds{cgr}{cs}(2:end);
    handles.cv=1;
    handles.cf=1;
    if ~isempty(get(handles.files_list,'String'))
        flist=flist(2:end);
        handles.dat.group(cgr).subject(cs).files=handles.dat.group(cgr).subject(cs).files(2:end,:);
    end
elseif cv==nvis
    vlist=vlist(1:end-1);
    handles.dat.group(cgr).subject(cs).visit=handles.dat.group(cgr).subject(cs).visit(1:end-1);
    handles.ds{cgr}{cs}=handles.ds{cgr}{cs}(1:end-1);
    handles.cv=handles.cv-1;
    handles.cf=1;
    if ~isempty(get(handles.files_list,'String'))
        flist=flist(1:end-1);
        handles.dat.group(cgr).subject(cs).files=handles.dat.group(cgr).subject(cs).files(1:end-1,:);
    end
else
    vlist=vlist([1:cv-1,cv+1:end]);
    handles.dat.group(cgr).subject(cs).visit=handles.dat.group(cgr).subject(cs).visit([1:cv-1,cv+1:end]);
    handles.ds{cgr}{cs}=handles.ds{cgr}{cs}([1:cv-1,cv+1:end]);
    handles.cv=handles.cv-1;
    handles.cf=1;
    if ~isempty(get(handles.files_list,'String'))
        flist=flist([1:cv-1,cv+1:end]);
        handles.dat.group(cgr).subject(cs).files=handles.dat.group(cgr).subject(cs).files([1:cv-1,cv+1:end],:);
    end
end
set(handles.visits_list,'String',vlist);
set(handles.files_list,'String',flist);
update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in file_add.
function file_add_Callback(hObject, eventdata, handles)
% hObject    handle to file_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%check that all previous fields were filled
if ~isfield(handles,'ds') || length(handles.ds)<1 || ...
        length(handles.ds{handles.cgr})<1 || ...
        length(handles.ds{handles.cgr}{handles.cs})<1
   beep
   disp('Please, select at least one group, subject and define all the visits before adding files')
   return
end
cgr=handles.cgr;
cs=handles.cs;
nv=length(handles.ds{handles.cgr}{handles.cs});
try
    prevlist=cellstr(handles.dat.group(cgr).subject(cs).files);
catch
    prevlist={};
end
fnames=spm_select(nv,'image','Select files for the subject in the same order used for visits',prevlist);
handles.dat.group(cgr).subject(cs).files=fnames;
handles.cf=1;
set(handles.files_list,'String',cellstr(fnames));
% Update handles structure
guidata(hObject, handles);

update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
% Update handles structure
guidata(hObject, handles);



% --- Executes on selection change in mask_menu.
function mask_menu_Callback(hObject, eventdata, handles)
% hObject    handle to mask_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mask_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mask_menu
val = get(hObject,'Value');
switch val
    case 1   % User wants no mask
        handles.dat.mask.val=0; % set flag "no mask"
        handles.dat.mask.fname={}; 
    case 2   % User wants a mask or the previous mask selected
        handles.dat.mask.val=1;
        tmp=get(handles.mask_menu,'String');
        if strcmpi(tmp{2},'yes')
            handles.dat.mask.fname=spm_select(1,'image','Select mask for the analysis');
            set(handles.mask_menu,'String',{'none';handles.dat.mask.fname; 'another'})
        else
            handles.dat.mask.fname=tmp{2};          
        end
    case 3  % User wants another mask
        handles.dat.mask.val=1;
        handles.dat.mask.fname=spm_select(1,'image','Select mask for the analysis');
        set(handles.mask_menu,'String',{'none';handles.dat.mask.fname; 'another'});
        set(hObject,'Value',2);
end    
handles.saved=0;
% Update handles structure
guidata(hObject, handles);
        

% --- Executes during object creation, after setting all properties.
function mask_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get and load swe.mat
swename = spm_select(1,'mat','Select swe.mat',[],pwd,'swe.mat');
swe=swe_load(swename);
if isempty(swe)
    beep
    disp('Could not load file')
    return
end
handles.dat=swe;
%Get the different fields and create the handles.ds cell array as well as
%complete the fields which might be missing (previous versions) and
%initialize the handles structure

handles.saved=1;
%flag if the masks are not linked to a modality name
if ~isfield(swe.masks,'mod_name')
    flagmask=0;
    handles.visits_list={};
else
    handles.visits_list={swe.masks(:).mod_name};
    flagmask=1;
end
    
%get the groups, subjects and modalities   
if isfield(swe,'group')
    ng=length(swe.group);
    handles.ds=cell(ng,1);
    for i=1:ng
        if isfield(swe.group(i),'subject')
            ns=length(swe.group(i).subject);
            handles.ds{i}=cell(1,ns);
            for j=1:ns
                if ~isfield(swe.group(i).subject,'subj_name')
                    handles.saved=0;
                    set(handles.save_data,'ForegroundColor',handles.color.high)
                    swe.group(i).subject(j).subj_name=['S',num2str(j)];
                end
                if isfield(swe.group(i).subject(j),'modality')
                    nm=length(swe.group(i).subject(j).modality);
                    handles.ds{i}{j}=cell(1,nm);
                    for k=1:nm
                        %if the flagmask is set to 0, then complete the list of
                        %modalities
                        if ~flagmask
                            mname=swe.group(i).subject(j).modality(k).mod_name;
                            if ~any(strcmpi(handles.visits_list,mname))
                                handles.saved=0;
                                handles.visits_list=[handles.visits_list, {mname}];
                            end
                        end
                        handles.ds{i}{j}{k}=size(swe.group(i).subject(j).modality(k).scans,1);
                    end
                end
            end
        end
    end
end

if ~flagmask==1
    disp('The files for masking are not linked to a modality name')
    disp('The information was erased. Please select the mask files')
    swe.masks=struct();
    handles.saved=0;
    set(handles.save_data,'ForegroundColor',handles.color.high)
end
if ~isempty(handles.visits_list)
    set(handles.mask_menu,'String',handles.visits_list);
else
    set(handles.mask_menu,'String',{'none'});
end
handles.cgr=1;
handles.cs=1;
handles.cv=1;
handles.cf=1;
a=fileparts(swename);
set(handles.edit1,'String',a)
set(handles.group_list,'String',{swe.group(:).gr_name})

%set the 'rename' and 'modify' right-clicks
%for groups
ren=uicontextmenu;
% Update handles structure
guidata(hObject, handles);
item1=uimenu(ren,'Label','Rename','Callback',@rengroup);
set(handles.group_list,'UIContextMenu',ren)
handles=guidata(hObject);
%for subjects
rens=uicontextmenu;
% Update handles structure
guidata(hObject, handles);
item1=uimenu(rens,'Label','Rename','Callback',@rensubj);
set(handles.subjects_list,'UIContextMenu',rens)
handles=guidata(hObject);
%for modalities
renm=uicontextmenu;
item1=uimenu(renm,'Label','Modify','Callback',@renmod);
% Update handles structure
guidata(hObject, handles);
set(handles.visits_list,'UIContextMenu',renm)
handles=guidata(hObject);

%Remove any field from previous computations
if isfield(swe,'fs')
    beep
    disp('Fields refering to feature sets have been found')
    disp('These will be removed if modifications to the dataset are performed')
    disp('Previously computed models will also be deleted')
    disp('Be sure to change the directory if you want to keep trace of previous work')
    handles.load_fsmod=1;
else
    handles.load_fsmod=0;
end

handles.dat=swe;
handles.dat.dir=a;
% Update handles structure
guidata(hObject, handles);

update_display_data(hObject,handles);
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in save_data.
function save_data_Callback(hObject, eventdata, handles)
% hObject    handle to save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.load_fsmod
    sure=swe_ui_sure;
    if ~sure || isempty(sure)
        return
    end
end

%Check that at least one group and one subject were entered
if ~isfield(handles,'ds') || length(handles.ds) <1 || ...
        length(handles.ds{1})<1 || length(handles.ds{1}{1})<1
    beep
    disp('Please, enter at least one group completed before saving')
    return
end

%Rearrange the data structure if the "Scans" option was selected for one of
%the groups
for i=1:length(handles.ds)
    if ~isempty(strfind(lower(handles.dat.group(i).subject(1).subj_name),'scan'))
        subj=struct();
        nsubj=handles.ds{i}{1}{1};
        for j=1:length(handles.dat.group(i).subject(1).modality)
            nsubj2=handles.ds{i}{1}{j};
            if nsubj ~=nsubj2
                beep
                sprintf('Number of subjects in modality %d and 1 of group %d are different ',j,i)
                disp('Please correct')
                return
            end
        end
        handles.ds{i}=cell(nsubj,1);
        for k=1:nsubj
            subj(k).subj_name=['S', num2str(k)];
            nmod=length(handles.dat.group(i).subject(1).modality);
            handles.ds{i}{k}=cell(nmod);
            for j=1:nmod
                handles.ds{i}{k}{j}=1;
                subj(k).modality(j)=handles.dat.group(i).subject(1).modality(j);
                subj(k).modality(j).scans=subj(k).modality(j).scans(k,:);
                if ~isempty(handles.dat.group(i).subject(1).modality(j).rt_subj)
                    subj(k).modality(j).rt_subj=subj(k).modality(j).rt_subj(k);
                else
                    subj(k).modality(j).rt_subj=[];
                end
                if ~isempty(handles.dat.group(i).subject(1).modality(j).covar)
                    subj(k).modality(j).covar=subj(k).modality(j).covar(k,:);
                else
                    subj(k).modality(j).covar=[];
                end
            end
        end
        handles.dat.group(i).subject=subj;
    end
end

%Check that the different groups have the same number of
%modalities and that the different subjects have the same number of
%modalities and files per modality
ng=length(handles.ds);
nm=length(handles.ds{1}{1});
ns=length(handles.ds{1});
list=get(handles.mask_menu,'String');
nmask=length(list);
%check that one mask was entered for each modality
if nmask~=length(handles.dat.masks)
    beep
    sprintf('%d masks were found, while %d modalities were added',length(handles.dat.masks),nmask)
    disp('Please correct')
    return
end

for i=1:ng
    matdat=zeros(ns,nm);
    ns=length(handles.ds{i});
    for j=1:ns
        nmi=length(handles.ds{i}{1});        
        nmj=length(handles.ds{i}{j});
        if nmj~=nmi
            beep
            sprintf('Numbers of modalities in subjects 1 and %d from group %d differ', j,i)
            disp('Please correct')
            return
        elseif nmj~=nm
            beep
            sprintf('Numbers of modalities in groups 1 and %d differ \n',i)
            disp('Please correct')
            return
        elseif nmj~=nmask
            beep
            sprintf('%d modalities found for subject %d of group %d, while %d masks found \n',nmj,j,i,nmask)
            disp('Possible errors in the modalities names, please correct')
            return
        end
        for k=1:nm
            m2=find(strcmpi({handles.dat.group(i).subject(j).modality(:).mod_name},list(k)));
            matdat(j,k)=handles.ds{i}{j}{m2};
            if isstruct(handles.dat.group(i).subject(j).modality(m2).design)
                des=handles.dat.group(i).subject(j).modality(m2).design;
                maxcond=max([des.conds(:).scans]);
                if matdat(j,k)<maxcond
                    beep
                    sprintf('Design of subject %d, group %d, modality %d, exceeds time series \n',j,i,k)
                    disp('Corresponding events were discarded')
                    for l=1:length(des.conds)
                        ovser=find(des.conds(l).scans>matdat(j,k));
                        inser=find(des.conds(l).scans<=matdat(j,k));
                        des.conds(l).discardedscans=[des.conds(l).discardedscans, des.conds(l).scans(ovser)];
                        des.conds(l).scans=des.conds(l).scans(inser);
                        des.conds(l).blocks=des.conds(l).blocks(inser);
                    end
                    handles.dat.group(i).subject(j).modality(m2).design=des;
                end
            end
        end
        handles.dat.group(i).subject(j).modality(k)=handles.dat.group(i).subject(j).modality(m2);
    end
end
def=swe_get_defaults('datad');
%save the data structure
disp('Saving the data.....>>')
swe=struct();
swe.group=handles.dat.group;

if ~isfield(handles.dat.group(1),'hrfoverlap')
    for i=1:ng
        swe.group(i).hrfoverlap=def.hrfw;
    end
end
if ~isfield(handles.dat.group(1),'hrfdelay')
    for i=1:ng
        swe.group(i).hrfdelay=def.hrfd;
    end
end
swe.masks=handles.dat.masks;

%Remove any field from previous computations if the swe is loaded and then
%modified
if isfield(swe,'fs')
    swe=rmfield(swe,'fs');
    beep
    disp('Fields refering to feature sets have been found')
    disp('These will be removed')
    disp('Be sure to change the directory if you want to keep trace of previous work')
end
if isfield(swe,'fas')
    swe=rmfield(swe,'fas');
end
if isfield(swe,'model')
    swe=rmfield(swe,'model');
end


resn=fullfile(handles.dat.dir,'swe.mat');
save(resn,'swe')

handles.saved=1;
disp('Save Done')
set(handles.save_data,'ForegroundColor',[0 0 0])

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in review_data.
function review_data_Callback(hObject, eventdata, handles)
% hObject    handle to review_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    nm=length(handles.dat.group(1).subject(1).modality);
catch
    beep
    disp('Please enter at least one subject in one group before reviewing')
    return
end
try
    nma=length(handles.dat.masks);
catch
    beep
    disp('A mask should be specified for each modality before reviewing')
    return
end
if nm~=nma
    beep
    disp('Number of masks does not match number of modalities')
    disp('Please, correct')
    return
end
if ~(handles.saved)
    swe=struct();
    swe.group=handles.dat.group;
    swe.masks=handles.dat.masks;
else
    fname=[get(handles.edit1,'String'),filesep,'swe.mat'];
    swe=swe_load(fname);
    if isempty(swe)
        beep
        disp('Could not load the saved swe.mat')
        return
    end
end
swe_data_review('UserData',{swe,handles.dat.dir});
fname=[get(handles.edit1,'String'),filesep,'swe.mat'];
swe=swe_load(fname);
handles.dat.group=swe.group;
handles.dat.mask=swe.mask;
% Update handles structure
guidata(hObject, handles);

           

% --- Executes on button press in quit_data.
function quit_data_Callback(hObject, eventdata, handles)
% hObject    handle to quit_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if no valuable information was entered, then exit
if ~isfield(handles,'ds') || length(handles.ds) <1 || ...
        length(handles.ds{1})<1 || length(handles.ds{1}{1})<1
    % The figure can be deleted now
    delete(handles.figure1);
else
    if ~handles.saved
        beep
        disp('Modifications performed since last saving, please save again')
        return
    else
        % The figure can be deleted now
        delete(handles.figure1);
    end
end

% --- Executes on selection change in swe_type_menu.
function swe_type_menu_Callback(hObject, eventdata, handles)
% hObject    handle to swe_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns swe_type_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from swe_type_menu
val = get(hObject,'Value');
list= get(hObject,'String');
handles.dat.swe_type=list{val};
handles.saved=0;
% Update handles structure
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function swe_type_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to swe_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in SS_menu.
function SS_menu_Callback(hObject, eventdata, handles)
% hObject    handle to SS_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SS_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SS_menu


% --- Executes during object creation, after setting all properties.
function SS_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SS_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
% --------------------- Subfunctions --------------------------------------
% -------------------------------------------------------------------------

function update_display_data(hObject,hand)

ng=hand.cgr;
ns=hand.cs;
nv=hand.cv;
nf=hand.cv;
%update group display
set(hand.group_list,'Value',ng)
%update subjects display
try
    set(hand.subjects_list,'String',{hand.dat.group(ng).subject(:).subj_name})
    set(hand.subjects_list,'Value',ns)
catch
    set(hand.subjects_list,'String',{})
end
%update visits display
try
    set(hand.visits_list,'String',{hand.dat.group(ng).subject(ns).visit(:).vis_name})
    set(hand.visits_list,'Value',nv)
catch
    set(hand.visits_list,'String',{})
end
%update file display
try  
    set(hand.files_list,'String',{hand.dat.group(ng).subject(ns).files})
    set(hand.files_list,'Value',nf)
catch
    set(hand.files_list,'String',{})
end
get(hand.files_list,'String')
% Update handles structure
guidata(hObject, hand);
