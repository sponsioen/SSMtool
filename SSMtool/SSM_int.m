function varargout = SSM_int(varargin)
% SSM_INT MATLAB code for SSM_int.fig
%      SSM_INT, by itself, creates a new SSM_INT or raises the existing
%      singleton*.
%
%      H = SSM_INT returns the handle to a new SSM_INT or the handle to
%      the existing singleton*.
%
%      SSM_INT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SSM_INT.M with the given input arguments.
%
%      SSM_INT('Property','Value',...) creates a new SSM_INT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SSM_int_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SSM_int_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SSM_int

% Last Modified by GUIDE v2.5 09-Mar-2017 08:06:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SSM_int_OpeningFcn, ...
                   'gui_OutputFcn',  @SSM_int_OutputFcn, ...
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

% --- Executes just before SSM_int is made visible.
function SSM_int_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SSM_int (see VARARGIN)

h = findobj('Tag','gui_main');

if ~isempty(h)
    handles.g1data = guidata(h);
end

resolution = get(0,'ScreenSize');
screen_width = resolution(3);
screen_height = resolution(4);
screen_ratio = screen_width/screen_height;
gui_int_ratio = 0.601713062098501;

if screen_ratio > 1
    screen_ref = screen_height;
else 
    screen_ref =  screen_width;
end

scale_par = 2.2;
gui_width = (screen_ref/scale_par)*gui_int_ratio;
gui_height = screen_ref/scale_par;
set(handles.gui_plot,'Units','pixels','Position',[handles.g1data.gui_main.Position(1),handles.g1data.gui_main.Position(2),gui_width,gui_height])
handles.spinner.position = [0,5,50,50];
% Choose default command line output for SSM_int
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = SSM_int_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function input_z20_Callback(hObject, eventdata, handles)
% hObject    handle to input_z20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: get(hObject,'String') returns contents of input_z20 as text
%        str2double(get(hObject,'String')) returns contents of input_z20 as a double


% --- Executes during object creation, after setting all properties.
function input_z20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_z20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function input_z10_Callback(hObject, eventdata, handles)
% hObject    handle to input_z10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: get(hObject,'String') returns contents of input_z10 as text
%        str2double(get(hObject,'String')) returns contents of input_z10 as a double


% --- Executes during object creation, after setting all properties.
function input_z10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_z10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function input_tend_Callback(hObject, eventdata, handles)
% hObject    handle to input_tend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_tend as text
%        str2double(get(hObject,'String')) returns contents of input_tend as a double


% --- Executes during object creation, after setting all properties.
function input_tend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_tend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in push_int.
function push_int_Callback(hObject, eventdata, handles)
% hObject    handle to push_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Spinner(handles.panel_spinner_int,'Busy','Start',handles.spinner.position);
mode = struct('WindowStyle','nonmodal','Interpreter','tex');
if (get(handles.checkbox_full_int,'Value')==1) && strcmp(get(handles.push_update,'Enable'),'on')
    txt = 'The parameterization coordinates have changed, please update your initial conditions for the full system.';
    h = errordlg(txt,'Input error',mode);
    Spinner(handles.panel_spinner_int,'Error','Stop',handles.spinner.position);
    return;
elseif (get(handles.checkbox_full_int,'Value')==1) && strcmp(get(handles.push_update,'Enable'),'off')
    try
        list_cor_Callback(hObject, eventdata, handles)
    catch
        Spinner(handles.panel_spinner_int,'Error','Stop',handles.spinner.position);
        return;
    end
end

h_gui = findobj('Tag','gui_plot');
handles = guidata(h_gui);
h = findobj('Tag','gui_main');

if ~isempty(h)
    handles.g1data = guidata(h);
end

if isempty(get(handles.input_z10,'String')) || isempty(get(handles.input_z20,'String')) || isempty(get(handles.input_tend,'String')) || isempty(get(handles.input_maxstep,'String')) 
    txt = 'Please fill in all the input fields.';
    h = errordlg(txt,'Input error',mode);
    Spinner(handles.panel_spinner_int,'Error','Stop',handles.spinner.position);
    return;
end

full_int = get(handles.checkbox_full_int,'Value');
z10 = str2double(get(handles.input_z10,'String'));
z20 = str2double(get(handles.input_z20,'String'));
tend = str2double(get(handles.input_tend,'String'));
options.maxstep = str2double(get(handles.input_maxstep,'String'));
options.complex_cor = handles.g1data.complex_cor;
y0 = [z10,z20];
sys.A = handles.g1data.sys.A;
sys.f = handles.g1data.sys.f;
sys.spv = handles.g1data.sys.spv;
y0_state_dum = handles.g1data.sys.Toutput;
nm = load('cs.mat');
SSM_function = @(z1,z2)eval(strcat('SSM_function_',nm.folder_id,'(z1,z2)'));

if full_int
    [t_full,ystate_full] = int_dyn(sys,tend,handles.y0_state_per,options);

    if handles.g1data.sys.modal
      ystate_full = (handles.g1data.sys.Vmodal\(ystate_full.')).';
    end
end

syms z1 z2
Zcor = [z1;z2];
[t,ystate] = int_red_dyn(tend,y0,Zcor,options);
R = handles.g1data.sys.Toutput;
Rho = ystate(:,1);
Theta = ystate(:,2);

if options.complex_cor
    [R{:}] = SSM_function(Rho.*exp(1i.*Theta),Rho.*exp(-1i.*Theta));
else
    [R{:}] = SSM_function(Rho.*cos(Theta),Rho.*sin(Theta));
end

if isfield(handles.g1data.plot,'g')
    delete(handles.g1data.plot.g)
end

if isfield(handles.g1data.plot,'w')
    delete(handles.g1data.plot.w)
end

handles.g1data.sys.R = R;
if ~isempty(handles.g1data.fig)
    if isvalid(handles.g1data.fig)
            figure(handles.g1data.fig)
            g = plot3(real(R{handles.g1data.plot.select(1)}),real(R{handles.g1data.plot.select(2)}),real(R{handles.g1data.plot.select(3)}));
            g.LineWidth = 2;
            g.Color = 'b';
            g.LineStyle = '-.';
            handles.g1data.plot.g = g;
            handles.g1data.legend_handle = {'SSM','Reduced Dynamics','Full System'};
            legend([handles.g1data.plot.h,g],handles.g1data.legend_handle{1:2})
    end
end
    
if full_int
    if ~isempty(handles.g1data.fig)
        if isvalid(handles.g1data.fig)
            figure(handles.g1data.fig)
            w = plot3(ystate_full(:,handles.g1data.plot.select(1)),ystate_full(:,handles.g1data.plot.select(2)),ystate_full(:,handles.g1data.plot.select(3)));
            w.LineWidth = 2.1;
            w.Color = 'r';
            w.LineStyle = '-';
            handles.g1data.plot.w = w;
            legend([handles.g1data.plot.h,g,w],handles.g1data.legend_handle)
        end
    end
end
Spinner(handles.panel_spinner_int,'Done','Stop',handles.spinner.position);
guidata(h,handles.g1data)
    
% --- Executes on button press in push_clear.
function push_clear_Callback(hObject, eventdata, handles)
% hObject    handle to push_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','gui_main');
Spinner(handles.panel_spinner_int,'Done','Stop',handles.spinner.position);
if ~isempty(h)
    handles.g1data = guidata(h);
end

if isfield(handles,'g1data')

    if isfield(handles.g1data,'legend_handle')
        handles.g1data.legend_handle =  {'SSM'};
        legend(handles.g1data.fig.CurrentAxes,handles.g1data.legend_handle)
    end

    if isfield(handles.g1data.plot,'g')
        delete(handles.g1data.plot.g)
    end

    if isfield(handles.g1data.plot,'w')
        delete(handles.g1data.plot.w)
    end
guidata(h,handles.g1data)
end
   
function input_maxstep_Callback(hObject, eventdata, handles)
% hObject    handle to input_maxstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_maxstep as text
%        str2double(get(hObject,'String')) returns contents of input_maxstep as a double


% --- Executes during object creation, after setting all properties.
function input_maxstep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_maxstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_init_pos(hObject, eventdata, handles)
mode = struct('WindowStyle','nonmodal','Interpreter','tex');
if isempty(get(handles.input_z10,'String')) || isempty(get(handles.input_z20,'String'))
    txt = 'Please fill in rho(0) and theta(0).';
    h = errordlg(txt,'Input error',mode);
    guidata(hObject, handles);
    Spinner(handles.panel_spinner_int,'Error','Stop',handles.spinner.position);
    error('Error');
end

z10 = str2double(get(handles.input_z10,'String'));
z20 = str2double(get(handles.input_z20,'String'));

if isnan(z10) || isnan(z20)
    txt = 'Please enter a numerical value for rho(0) and theta(0).';
    h = errordlg(txt,'Input error',mode);
    guidata(hObject, handles);
    Spinner(handles.panel_spinner_int,'Error','Stop',handles.spinner.position);
    error('Error');
end
 
set(handles.list_cor,'Enable','On');
set(handles.list_cor,'String',{handles.g1data.sys.Toutput{:}});
set(handles.push_update,'Enable','Off');
nm = load('cs.mat');
SSM_function = @(z1,z2)eval(strcat('SSM_function_',nm.folder_id,'(z1,z2)'));
y0 = [z10,z20];
y0_state_dum = handles.g1data.sys.Toutput;
    
if handles.g1data.sys.modal
    if handles.g1data.complex_cor
        [y0_state_dum{:}] = SSM_function(y0(1)*exp(1i*y0(2)),y0(1)*exp(-1i.*y0(2)));
    else
        [y0_state_dum{:}] = SSM_function(y0(1)*cos(y0(2)),y0(1)*sin(y0(2)));
    end
    y0_state = real(handles.g1data.sys.Vmodal*([y0_state_dum{:}].')).';
    handles.y0_state_per_modal = real([y0_state_dum{:}]);
    handles.y0_state_per = y0_state;
    set(handles.input_per,'Enable','On');
    set(handles.input_per,'String',num2str(handles.y0_state_per_modal(get(handles.list_cor,'Value')),'%.14f'));
else
    if handles.g1data.complex_cor
        [y0_state_dum{:}] = SSM_function(y0(1)*exp(1i*y0(2)),y0(1)*exp(-1i.*y0(2)));
    else
        [y0_state_dum{:}] = SSM_function(y0(1)*cos(y0(2)),y0(1)*sin(y0(2)));
    end
    y0_state = real([y0_state_dum{:}]);
    handles.y0_state_per = y0_state;
    set(handles.input_per,'Enable','On');
    set(handles.input_per,'String',num2str(y0_state(get(handles.list_cor,'Value')),'%.14f'));
end
    
    guidata(hObject, handles);

% --- Executes on button press in checkbox_full_int.
function checkbox_full_int_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_full_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    try
        update_init_pos(hObject, eventdata, handles)
    catch 
        Spinner(handles.panel_spinner_int,'Error','Stop',handles.spinner.position);
        set(handles.checkbox_full_int,'Value',0);
    end
else
    set(handles.list_cor,'Enable','Off');
    set(handles.push_update,'Enable','Off');
    set(handles.input_per,'Enable','Off');
end

% --- Executes on selection change in list_cor.
function list_cor_Callback(hObject, eventdata, handles)
% hObject    handle to list_cor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

oldValues = get(handles.list_cor,'UserData');
newValues = get(handles.list_cor,'Value');
if numel(newValues) > 1    
    newValues = oldValues;
end  
set(handles.list_cor,'Value',newValues)
set(handles.list_cor,'UserData',newValues)
mode = struct('WindowStyle','nonmodal','Interpreter','tex');
if isempty(get(handles.input_per,'String')) 
    txt = 'Please fill in an initial coordinate.';
    h = errordlg(txt,'Input error',mode);
    guidata(hObject, handles);
    set(handles.list_cor,'Value',oldValues);
    set(handles.list_cor,'UserData',oldValues);
    Spinner(handles.panel_spinner_int,'Error','Stop',handles.spinner.position);
    error('Please fill in an initial coordinate.');
    
end
    
if isnan(str2double(get(handles.input_per,'String'))) 
    txt = 'Please enter a numerical value as an initial coordinate.';
    h = errordlg(txt,'Input error',mode);
    guidata(hObject, handles);
    set(handles.list_cor,'Value',oldValues);
    set(handles.list_cor,'UserData',oldValues);
    Spinner(handles.panel_spinner_int,'Error','Stop',handles.spinner.position);
    error('Please enter a numerical value as an initial coordinate.');
end
    
if handles.g1data.sys.modal
    handles.y0_state_per_modal(oldValues) = str2double(get(handles.input_per,'String'));
    handles.y0_state_per = real(handles.g1data.sys.Vmodal*(handles.y0_state_per_modal.')).';
    set(handles.input_per,'String',num2str(handles.y0_state_per_modal(get(handles.list_cor,'Value')),'%.14f'));
else 
    handles.y0_state_per(oldValues) = str2double(get(handles.input_per,'String'));
    set(handles.input_per,'String',num2str(handles.y0_state_per(get(handles.list_cor,'Value')),'%.14f'));
end
   
guidata(hObject, handles);
  
% Hints: contents = cellstr(get(hObject,'String')) returns list_cor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_cor


% --- Executes during object creation, after setting all properties.
function list_cor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_cor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function input_per_Callback(hObject, eventdata, handles)
% hObject    handle to input_per (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_per as text
%        str2double(get(hObject,'String')) returns contents of input_per as a double


% --- Executes during object creation, after setting all properties.
function input_per_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_per (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on key press with focus on input_z10 and none of its controls.
function input_z10_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to input_z10 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if get(handles.checkbox_full_int,'Value')
    set(handles.push_update,'Enable','On');
    set(handles.list_cor,'Enable','Off');
    set(handles.input_per,'Enable','Off');
end

% --- Executes on button press in push_update.
function push_update_Callback(hObject, eventdata, handles)
% hObject    handle to push_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 try
    update_init_pos(hObject, eventdata, handles)
 catch ME
     Spinner(handles.panel_spinner_int,'Error','Stop',handles.spinner.position);
     rethrow(ME);
 end

% --- Executes on key press with focus on input_z20 and none of its controls.
function input_z20_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to input_z20 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox_full_int,'Value')
    set(handles.push_update,'Enable','On');
    set(handles.list_cor,'Enable','Off');
    set(handles.input_per,'Enable','Off');
end
