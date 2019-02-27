function varargout = SSM_exp(varargin)
% SSM_EXP MATLAB code for SSM_exp.fig
%      SSM_EXP, by itself, creates a new SSM_EXP or raises the existing
%      singleton*.
%
%      H = SSM_EXP returns the handle to a new SSM_EXP or the handle to
%      the existing singleton*.
%
%      SSM_EXP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SSM_EXP.M with the given input arguments.
%
%      SSM_EXP('Property','Value',...) creates a new SSM_EXP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SSM_exp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SSM_exp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SSM_exp

% Last Modified by GUIDE v2.5 01-May-2017 11:33:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SSM_exp_OpeningFcn, ...
                   'gui_OutputFcn',  @SSM_exp_OutputFcn, ...
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


% --- Executes just before SSM_exp is made visible.
function SSM_exp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SSM_exp (see VARARGIN)

handles.fig_backbone = [];
handles.plot = [];
h = findobj('Tag','gui_main');

if ~isempty(h)
handles.g1data = guidata(h);
end

resolution = get(0,'ScreenSize');
screen_width = resolution(3);
screen_height = resolution(4);
screen_ratio = screen_width/screen_height;
gui_exp_ratio = 1.2168;

if screen_ratio > 1
    screen_ref = screen_height;
else 
    screen_ref =  screen_width;
end
    
scale_par = 1.4;
gui_width = (screen_ref/scale_par)*gui_exp_ratio;
gui_height = screen_ref/scale_par;
set(handles.gui_exp,'Units','pixels','Position',[handles.g1data.gui_main.Position(1),handles.g1data.gui_main.Position(2),gui_width,gui_height]); 
nm = load('cs.mat');
SSM_function = @(z1,z2)eval(strcat('SSM_function_',nm.folder_id,'(z1,z2)'));
R_function = @(z1,z2)eval(strcat('R_function_',nm.folder_id,'(z1,z2)'));
n_dof = numel(handles.g1data.sys.spv);
modal =  handles.g1data.sys.modal;
complex = handles.g1data.sys.complex;
T = cell(1,n_dof);   

if modal 
    for i = 1:n_dof
       T{i} = sprintf('q_%d',i);
    end
elseif complex 
    for i = 1:n_dof
       T{i} = sprintf('c_%d',i);
    end
else 
    for i = 1:n_dof
       T{i} = char(handles.g1data.sys.spv(i));
    end
end
    
X_cor = cell(1,n_dof/2); 
for i = 1:n_dof/2
   X_cor{i} = strcat(strcat('|',char(handles.g1data.sys.spv(i))),'|');
end

set(handles.popup_cor_avg,'String',{X_cor{:}});
set(handles.popup_cor_avg,'Value',1);
handles.spinner.position = [0,5,50,50];
syms z1 z2
S = T;

if handles.g1data.sys.int_res || handles.g1data.sys.ho_int_res
    set(handles.input_par,'Enable','On')
    set(handles.input_points,'Enable','On')
    set(handles.push_backbone,'Enable','On')
    set(handles.checkbox_avg,'Enable','On');
    set(handles.popup_cor_avg,'Enable','On');
else
    set(handles.input_par,'Enable','Off')
    set(handles.input_points,'Enable','Off')
    set(handles.push_backbone,'Enable','Off')
    set(handles.checkbox_avg,'Enable','Off');
    set(handles.popup_cor_avg,'Enable','Off');
end

[S{:}] = SSM_function(z1,z2);
R = {'R1','R2'};
Rpol = {'r','theta'};
syms  r theta real

if handles.g1data.complex_cor
    [R{:}] = R_function(r*exp(1i*theta),r*exp(-1i*theta));
    alpha = real(simplify(R{1}));
    beta = imag(simplify(R{1}));
else
    [R{:}] = R_function(r*cos(theta),r*sin(theta));
    alpha = real(R{1});
    beta = real(R{2});
end

R{1} = expand(simplify(alpha*cos(theta) + beta*sin(theta)));
R{2} = expand(simplify((1/r)*(beta*cos(theta)-alpha*sin(theta))));   % + or - 

if numel(symvar(jacobian(R{2},theta)))>1
    handles.back_avg = 0;
    handles.mono_var_freq  = 0;

    set(handles.checkbox_avg,'Value',0);
    set(handles.checkbox_avg,'Enable','Off');
    set(handles.popup_cor_avg,'Enable','On');
else   
    handles.mono_var_freq  = 1;
end  

SSM_expression = cell(2*n_dof,1);
R_expression = cell(4,1);
   
m = 1;
for j = 1:2:2*n_dof
    SSM_expression{j} = sprintf('%s(z1,z2) = %s',char(T{m}),char(vpa(S{m},5)));
    SSM_expression{j+1} = char('');
    m = m+1;
end

n = 1;
for j = 1:2:4
    R_expression{j} = sprintf('d%s/dt = %s',Rpol{n},char(vpa(R{n},5)));
    R_expression{j+1} = char('');
    n= n+1;
end
   
set(handles.text_SSM,'String',SSM_expression)
set(handles.text_R,'String',R_expression)

% Choose default command line output for SSM_exp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SSM_exp wait for user response (see UIRESUME)
% uiwait(handles.gui_exp);


% --- Outputs from this function are returned to the command line.
function varargout = SSM_exp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function text_SSM_Callback(hObject, eventdata, handles)
% hObject    handle to text_SSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_SSM as text
%        str2double(get(hObject,'String')) returns contents of text_SSM as a double


% --- Executes on button press in push_backbone.
function push_backbone_Callback(hObject, eventdata, handles)
% hObject    handle to push_backbone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Spinner(handles.panel_spinner_exp,'Busy','Start',handles.spinner.position);
syms rho theta real
nm = load('cs.mat');
SSM_function = @(z1,z2)eval(strcat('SSM_function_',nm.folder_id,'(z1,z2)'));
input_par = str2num(get(handles.input_par,'String'));
input_points = str2num(get(handles.input_points,'String'));
r_min = 0.00001;
mode = struct('WindowStyle','nonmodal','Interpreter','tex');
if isempty(input_par) || isempty(input_points)
    txt = 'Please fill in all required fields';
    h = errordlg(txt,'Input error',mode);
    Spinner(handles.panel_spinner_exp,'Error','Stop',handles.spinner.position);
    return;
elseif input_par <= r_min
    txt = 'Parameterization domain is too small';
    h = errordlg(txt,'Input error',mode);
    Spinner(handles.panel_spinner_exp,'Error','Stop',handles.spinner.position);
    return;
end

S = cell(1,numel(handles.g1data.sys.spv));

if handles.g1data.complex_cor
    [S{:}] = SSM_function(rho.*exp(1i.*theta),rho.*exp(-1i.*theta));
else
    [S{:}] = SSM_function(rho.*cos(theta),rho.*sin(theta));
end

if handles.g1data.sys.modal
    S_vec = handles.g1data.sys.Vmodal*([S{:}].');
elseif handles.g1data.sys.complex
    S_vec = handles.g1data.sys.V*([S{:}].');
else
    S_vec = [S{:}].';
end

handles.back_avg = get(handles.checkbox_avg,'Value');


if handles.back_avg
    handles.mono_var_freq = 1;
    h = @(rho,theta)real(eval(sqrt(S_vec(1:2).'*S_vec(1:2))));
    amp = @(rho)(1/(2*pi)).*integral(@(theta)h(rho,theta),0,2*pi);
else
    val = get(handles.popup_cor_avg,'Value');
    h = @(rho,theta)real(eval(S_vec(val)));
    theta_linspace = linspace(0,2*pi,200);
    amp_vec = @(rho) abs(h(rho,theta_linspace));
end
R = cell(1,2);
if handles.g1data.complex_cor
    [R{:}] = R_sub_function(rho*exp(1i*theta),rho.*exp(-1i.*theta));
    alpha = real(R{1});
    beta = imag(R{1});
else
    [R{:}] = R_sub_function(rho.*cos(theta),rho.*sin(theta));
    alpha = real(R{1});
    beta = real(R{2});
end

R{1} = simplify(alpha.*cos(theta) + beta.*sin(theta));
R{2} = simplify((1/rho).*(beta.*cos(theta)-alpha.*sin(theta)));

if handles.mono_var_freq    
    omega = @(rho)eval(R{2});
else
    omega = @(rho,theta)eval(R{2});
end

rho_dummy = linspace(r_min,input_par,input_points);
amp_plot = zeros(1,numel(rho_dummy));
theta_plot = zeros(1,numel(rho_dummy));

if handles.mono_var_freq  
    if handles.back_avg
        for i=1:numel(rho_dummy)
            amp_plot(i) =  amp(rho_dummy(i));
        end
    else
        for i=1:numel(rho_dummy)
            amp_plot(i) = max(amp_vec(rho_dummy(i)));
        end
    end
else
    for i=1:numel(rho_dummy)
        [amp_max,theta_i] = max(amp_vec(rho_dummy(i)));
        theta_plot(i) = theta_linspace(theta_i);
        amp_plot(i) =  amp_max;
    end
end

if isfield(handles.plot,'s')
    delete(handles.plot.s)
end

if ~isempty(handles.fig_backbone)
    if isvalid(handles.fig_backbone)
        clf(handles.fig_backbone)
        close(handles.fig_backbone)
    end
end   
            
handles.fig_backbone = figure;
if ~handles.mono_var_freq  
    s = plot(omega(rho_dummy,theta_plot),amp_plot,'-','LineWidth',2);
else
    s = plot(omega(rho_dummy),amp_plot,'-','LineWidth',2);
end
handles.axis_handle_backbone = handles.fig_backbone.CurrentAxes;
handles.plot.s = s;
set(handles.axis_handle_backbone,'FontSize',15)
xlabel('Frequency (rad/s)')
ylabel('Amplitude')
grid on
axis auto
axis square
box on
set(handles.axis_handle_backbone,'FontUnits','normalized');
Spinner(handles.panel_spinner_exp,'Done','Stop',handles.spinner.position);

% Update handles structure
guidata(hObject, handles);



function input_par_Callback(hObject, eventdata, handles)
% hObject    handle to input_par (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_par as text
%        str2double(get(hObject,'String')) returns contents of input_par as a double


% --- Executes during object creation, after setting all properties.
function input_par_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_par (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_points_Callback(hObject, eventdata, handles)
% hObject    handle to input_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_points as text
%        str2double(get(hObject,'String')) returns contents of input_points as a double


% --- Executes during object creation, after setting all properties.
function input_points_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_avg.
function checkbox_avg_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_avg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') 
   set(handles.popup_cor_avg,'Enable','Off')
   handles.back_avg = 0;
else
   set(handles.popup_cor_avg,'Enable','On')
   handles.back_avg = 1;
end

% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox_avg


% --- Executes on selection change in popup_cor_avg.
function popup_cor_avg_Callback(hObject, eventdata, handles)
% hObject    handle to popup_cor_avg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_cor_avg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_cor_avg


% --- Executes during object creation, after setting all properties.
function popup_cor_avg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_cor_avg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
