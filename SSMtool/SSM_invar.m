function varargout = SSM_invar(varargin)
% SSM_INVAR MATLAB code for SSM_invar.fig
%      SSM_INVAR, by itself, creates a new SSM_INVAR or raises the existing
%      singleton*.
%
%      H = SSM_INVAR returns the handle to a new SSM_INVAR or the handle to
%      the existing singleton*.
%
%      SSM_INVAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SSM_INVAR.M with the given input arguments.
%
%      SSM_INVAR('Property','Value',...) creates a new SSM_INVAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SSM_invar_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SSM_invar_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SSM_invar

% Last Modified by GUIDE v2.5 31-Aug-2017 17:34:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SSM_invar_OpeningFcn, ...
                   'gui_OutputFcn',  @SSM_invar_OutputFcn, ...
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


% --- Executes just before SSM_invar is made visible.
function SSM_invar_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SSM_invar (see VARARGIN)

h = findobj('Tag','gui_main');

if ~isempty(h)
    handles.g1data = guidata(h);
end

handles.spinner.position = [0,5,50,50];
resolution = get(0,'ScreenSize');
screen_width = resolution(3);
screen_height = resolution(4);
screen_ratio = screen_width/screen_height;
gui_int_ratio = 0.83938224;

if screen_ratio > 1
    screen_ref = screen_height;
else 
    screen_ref =  screen_width;
end

scale_par = 3.2;
gui_width = (screen_ref/scale_par)*gui_int_ratio;
gui_height = screen_ref/scale_par;
set(handles.gui_invar,'Units','pixels','Position',[handles.g1data.gui_main.Position(1),handles.g1data.gui_main.Position(2),gui_width,gui_height])

% Choose default command line output for SSM_invar
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SSM_invar wait for user response (see UIRESUME)
% uiwait(handles.gui_invar);


% --- Outputs from this function are returned to the command line.
function varargout = SSM_invar_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_rho_0_Callback(hObject, eventdata, handles)
% hObject    handle to input_rho_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_rho_0 as text
%        str2double(get(hObject,'String')) returns contents of input_rho_0 as a double


% --- Executes during object creation, after setting all properties.
function input_rho_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_rho_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_rho_1_Callback(hObject, eventdata, handles)
% hObject    handle to input_rho_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_rho_1 as text
%        str2double(get(hObject,'String')) returns contents of input_rho_1 as a double


% --- Executes during object creation, after setting all properties.
function input_rho_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_rho_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_N_Callback(hObject, eventdata, handles)
% hObject    handle to input_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_N as text
%        str2double(get(hObject,'String')) returns contents of input_N as a double


% --- Executes during object creation, after setting all properties.
function input_N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_T_Callback(hObject, eventdata, handles)
% hObject    handle to input_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_T as text
%        str2double(get(hObject,'String')) returns contents of input_T as a double


% --- Executes during object creation, after setting all properties.
function input_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_start.
function push_start_Callback(hObject, eventdata, handles)
% hObject    handle to push_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_result,'String','');
Spinner(handles.panel_spinner_invar,'Busy','Start',handles.spinner.position);
if isempty(get(handles.input_rho_0,'String')) || isempty(get(handles.input_rho_1,'String')) || isempty(get(handles.input_N,'String')) || isempty(get(handles.input_T,'String')) 
    txt = 'Please fill in all the input fields.';
    h = errordlg(txt,'Input error',mode);
    return;
end

rho_0 = str2double(get(handles.input_rho_0,'String'));
rho_1 = str2double(get(handles.input_rho_1,'String'));
N = str2double(get(handles.input_N,'String'));
T = str2double(get(handles.input_T,'String'));

try
    error = measure_inv_autonomous(rho_0,rho_1,N,T);
    set(handles.text_result,'String',sprintf('%d',error));
catch ME
    Spinner(handles.panel_spinner_invar,'Error','Stop',handles.spinner.position);
    rethrow(ME);
end
Spinner(handles.panel_spinner_invar,'Done','Stop',handles.spinner.position);



function text_result_Callback(hObject, eventdata, handles)
% hObject    handle to text_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_result as text
%        str2double(get(hObject,'String')) returns contents of text_result as a double


% --- Executes during object creation, after setting all properties.
function text_result_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
