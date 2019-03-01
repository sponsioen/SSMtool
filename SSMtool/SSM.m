function varargout = SSM(varargin)
% SSM MATLAB code for SSM.fig
%      SSM, by itself, creates a new SSM or raises the existing
%      singleton*.
%
%      H = SSM returns the handle to a new SSM or the handle to
%      the existing singleton*.
%
%      SSM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SSM.M with the given input arguments.
%
%      SSM('Property','Value',...) creates a new SSM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SSM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SSM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SSM

% Last Modified by GUIDE v2.5 01-Sep-2017 10:04:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SSM_OpeningFcn, ...
                   'gui_OutputFcn',  @SSM_OutputFcn, ...
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


% --- Executes just before SSM is made visible.
function SSM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SSM (see VARARGIN)
    
    clc
    warning('off', 'symbolic:sym:sym:DeprecateExpressions')
    
    if exist(strcat(pwd,'/Data'),'dir')==7
        addpath(strcat(pwd,'/Data/'));
        rmpath(genpath(strcat(pwd,'/Data/')));
        addpath(genpath(strcat(pwd,'/Data/')));
    else
        mkdir(strcat(pwd,'/Data'))
        addpath(strcat(pwd,'/Data'))
    end
    
    handles.plot.xmin = str2num(get(handles.plot_xmin,'String'));
    handles.plot.xmax = str2num(get(handles.plot_xmax,'String'));

    handles.plot.ymin = str2num(get(handles.plot_ymin,'String'));
    handles.plot.ymax = str2num(get(handles.plot_ymax,'String'));

    handles.plot.zmin = str2num(get(handles.plot_zmin,'String'));
    handles.plot.zmax = str2num(get(handles.plot_zmax,'String'));

    handles.plot.r_max = str2num(get(handles.plot_r_max,'String'));
    handles.fig = [];
    handles.plot.cor = 1;
    handles.plot.select = 1:3;
    handles.spinner.position = [0,5,50,50];
    
    resolution = get(0,'ScreenSize');
    screen_width = resolution(3);
    screen_height = resolution(4);
    screen_ratio = screen_width/screen_height;
    gui_ratio = 1.1925;
    
    if screen_ratio > 1
        screen_ref = screen_height;
    else 
        screen_ref =  screen_width;
    end
    
    scale_par = 1.2;
    gui_width = (screen_ref/scale_par)*gui_ratio;
    gui_height = screen_ref/scale_par;

    set(handles.gui_main,'Units','pixels','Position',[0,0,gui_width,gui_height])
    
    
    disp('----------------------------------------------------------------')
    disp('Automated computation of Spectral Submanifolds for autonomous')
    disp('mechanical systems V1.0 (2017)')
    disp(' ')
    disp('Sten Ponsioen (stenp@ethz.ch)')
    disp('George Haller (georgehaller@ethz.ch)')
    disp('----------------------------------------------------------------')

% Choose default command line output for SSM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SSM wait for user response (see UIRESUME)
% uiwait(handles.gui_main);


% --- Outputs from this function are returned to the command line.
function varargout = SSM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popup_academic.
function popup_academic_Callback(hObject, eventdata, handles)
% hObject    handle to popup_academic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_academic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_academic


% --- Executes during object creation, after setting all properties.
function popup_academic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_academic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_order.
function popup_order_Callback(hObject, eventdata, handles)
% hObject    handle to popup_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: contents = cellstr(get(hObject,'String')) returns popup_order contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_order


% --- Executes during object creation, after setting all properties.
function popup_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_load.
function push_load_Callback(hObject, eventdata, handles)
% hObject    handle to push_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function input_mass_Callback(hObject, eventdata, handles)
% hObject    handle to input_mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_mass as text
%        str2double(get(hObject,'String')) returns contents of input_mass as a double


% --- Executes during object creation, after setting all properties.
function input_mass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_nonlinear_Callback(hObject, eventdata, handles)
% hObject    handle to input_nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_nonlinear as text
%        str2double(get(hObject,'String')) returns contents of input_nonlinear as a double


% --- Executes during object creation, after setting all properties.
function input_nonlinear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_nonlinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_damping_Callback(hObject, eventdata, handles)
% hObject    handle to input_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_damping as text
%        str2double(get(hObject,'String')) returns contents of input_damping as a double


% --- Executes during object creation, after setting all properties.
function input_damping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_stiffness_Callback(hObject, eventdata, handles)
% hObject    handle to input_stiffness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_stiffness as text
%        str2double(get(hObject,'String')) returns contents of input_stiffness as a double


% --- Executes during object creation, after setting all properties.
function input_stiffness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_stiffness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_pos_Callback(hObject, eventdata, handles)
% hObject    handle to input_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_pos as text
%        str2double(get(hObject,'String')) returns contents of input_pos as a double


% --- Executes during object creation, after setting all properties.
function input_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_vel_Callback(hObject, eventdata, handles)
% hObject    handle to input_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_vel as text
%        str2double(get(hObject,'String')) returns contents of input_vel as a double


% --- Executes during object creation, after setting all properties.
function input_vel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_load_ms.
function push_load_ms_Callback(hObject, eventdata, handles)
% hObject    handle to push_load_ms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    clear_spectral_handles(hObject, eventdata, handles)
  
    if isfield(handles,'sys')
        handles = rmfield(handles,'sys');
    end

    Spinner(handles.panel_spinner,'Busy','Start',handles.spinner.position); 
    handles.sys.symbolic = 0;
    handles.sys.scaling = str2double(get(handles.scaling_factor,'String'));
    radio_select_con = get(handles.buttongroup_con,'SelectedObject');
    
    switch radio_select_con.String
        case 'Yes'
             handles.sys.conservative = 1;
        case 'No'
             handles.sys.conservative = 0;
    end
 
    M_string = get(handles.input_mass,'String'); 
    C_string = get(handles.input_damping,'String'); 
    K_string = get(handles.input_stiffness,'String'); 
    
    q_state = get(handles.input_pos,'String');
    qd_state = get(handles.input_vel,'String');
    
    q_state(q_state==' ') = '';
    q_state([strfind(q_state,'['),strfind(q_state,']')]) = [];
    q = sym(strread(q_state,'%s','delimiter',';'));

    qd_state(qd_state==' ') = '';
    qd_state([strfind(qd_state,'['),strfind(qd_state,']')]) = [];
    qd = sym(strread(qd_state,'%s','delimiter',';'));
    
    if isempty(q_state) || isempty(qd_state)
        errordlg('State variables are not specified','Input Error');
        Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
    elseif isempty(M_string) || isempty(C_string) || isempty(K_string)
        errordlg('System matrices are not specified','Input Error');
        Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
    elseif isempty(get(handles.input_nonlinear,'String'))
        errordlg('Nonlinear vector is not specified','Input Error');
        Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
    elseif isempty(get(handles.scaling_factor,'String'))
        errordlg('Scaling factor is not specified','Input Error');
        Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
    else
        
    set(handles.push_clear_loaded,'Enable','On');

    spv  = [q;qd];
    
    for i=1:numel(spv)
       eval(strcat(char(spv(i)),sprintf('= spv(%d);',i)));
    end

    f = eval(get(handles.input_nonlinear,'String'));
    
    MV = ver('MATLAB');
    
    if sum(MV.Version == '9.1') == 3
            M =  sym(M_string);
    else
            M =  str2sym(M_string);
    end
    
    if handles.sys.conservative
        C = zeros(numel(q));
    else
        if sum(MV.Version == '9.1') == 3
            C = sym(C_string);
        else   
            C = str2sym(C_string);   
        end
    end
    
    if sum(MV.Version == '9.1') == 3
        K =  sym(K_string);
    else
        K =  str2sym(K_string);
    end
    
    fnl = [sym(zeros(numel(q),1));-M\f];
    
    try
        [lambda,T,A] = compute_subspace(M,C,K,handles.sys.scaling,handles.sys.conservative);
    catch NE
        Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
        rethrow(NE)
    end

    handles.sys.M = M;
    handles.sys.C = C;
    handles.sys.K = K;
    handles.sys.lambda = lambda;
    handles.sys.T = T;
    handles.sys.A = A;
    handles.sys.f = fnl;
    handles.sys.spv = spv;

    set(handles.list_lambda,'String',num2str(lambda.num));    
    set(handles.push_select_lambda,'Enable','On');
    
    handles.compute.lambda = 0;
    handles.compute.order = 1;
    set(handles.text_order,'String','5');
    set(handles.text_order,'Enable','Off');
   
    check_status_compute(hObject, eventdata, handles)
    set(handles.text_file_loaded,'String', '')
    set(handles.text_file_loaded,'Visible','Off')
    Spinner(handles.panel_spinner,'Done','Stop',handles.spinner.position);
    
    % Update handles structure
    guidata(hObject, handles);
    end
catch ME
    errordlg('Input error mechanical system.','Input Error');
    Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
    rethrow(ME);
end


% --- Executes on selection change in list_lambda.
function list_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to list_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    oldValues = get(handles.list_lambda,'UserData');
    newValues = get(handles.list_lambda,'Value');
    if numel(newValues) > 2    
        newValues = oldValues;
    end
    set(handles.list_lambda,'Value',newValues)
    set(handles.list_lambda,'UserData',newValues)
    
  % Update handles structure
  guidata(hObject, handles);
    

% --- Executes during object creation, after setting all properties.
function list_lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_select_lambda.
function push_select_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to push_select_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.text_nonres,'String','');
    set(handles.text_nonres_2,'String','');
    set(handles.text_sigma,'String','');

    getValues = get(handles.list_lambda,'Value');
 
        if numel(getValues) < 2    
            mode = struct('WindowStyle','nonmodal','Interpreter','tex');
            h = errordlg('Please select two eigenvalues.','Selection Error', mode);
            handles.compute.lambda = 0;
            check_status_compute(hObject, eventdata, handles)      
        else
            Spinner(handles.panel_spinner,'Busy','Start',handles.spinner.position);
            
            A = handles.sys.A;
            lambda =  handles.sys.lambda.sym;
            diff = 1:numel(lambda); 
            p = setdiff(diff,getValues);
 
            handles.sys.lambda_select = lambda(getValues);
            lambda_remain = lambda(p);
            handles.sys.lambda_remain = lambda_remain;
            handles.sys.lambda_p = p;
            handles.sys.lambda_s = getValues;
            
            try
                handles_out = orderT(A,handles);
            catch Me
                rethrow(Me)
            end

            try
                handles = check_res(handles_out); 
            catch Me
                handles.compute.lambda = 0;
                check_status_compute(hObject, eventdata, handles)
                rethrow(Me)
            end
                    
            handles.compute.lambda = 1;
            check_status_compute(hObject, eventdata, handles)
            
            set(handles.text_order,'Enable','On');
            
            if ~handles.sys.conservative
                set(handles.check_higher_int_res,'Enable','On');
            end

            % Update handles structure
            guidata(hObject, handles);
            Spinner(handles.panel_spinner,'Done','Stop',handles.spinner.position);
            
        end
        
        
function check_status_compute(hObject, eventdata, handles)

    handles.compute.order = ~isnan(str2double(get(handles.text_order,'String')));
 
    if handles.compute.lambda && handles.compute.order
       set(handles.push_compute,'Enable','On')
    else
       set(handles.push_compute,'Enable','Off')
    end
       
    % Update handles structure
    guidata(hObject, handles);
        
        
        
        

% --- Executes on button press in push_compute.
function push_compute_Callback(hObject, eventdata, handles)
% hObject    handle to push_compute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    Spinner(handles.panel_spinner,'Busy','Start',handles.spinner.position);
    radio_select = get(handles.buttongroup_cor,'SelectedObject');
    mode = struct('WindowStyle','nonmodal','Interpreter','tex');
    handles.sys.n = str2double(get(handles.text_order,'String'));
    
    if isnan(handles.sys.n)     
        txt = 'Input error for order of SSM expansion';
        h = errordlg(txt,'Input Error',mode);
        Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
        return;
    elseif handles.sys.n>50
        txt = 'The maximum order of SSM expansion is currently set to 50';
        h = errordlg(txt,'Order Error',mode);
        Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
        return;
    end   
    
    c = fix(clock);
    folder_id = sprintf('%d_%d_%d_%dh_%dm_%ds',c([3,2,1,4,5,6]));
    save('cs.mat','folder_id');
    
    
    switch radio_select.String
        case 'Physical'
             handles.sys.modal = 0;
             handles.sys.complex = 0;
        case 'Modal'
             handles.sys.modal = 1;
             handles.sys.complex = 0;
        case 'Complex'
             handles.sys.modal = 0;
             handles.sys.complex = 1;
    end

    if ~isempty(handles.fig)
        if isvalid(handles.fig)
            clf(handles.fig)
            close(handles.fig)
        end
    end
    

    h1 = findobj('Tag','gui_plot');
    h2 = findobj('Tag','gui_exp');
    
    if ~isempty(h1)
        handles.gui_plot_data = guidata(h1);
        close(handles.gui_plot_data.gui_plot)
    end
    
    if ~isempty(h2)
        handles.gui_exp_data = guidata(h2);
        
         if ~isempty(handles.gui_exp_data.fig_backbone)
            if isvalid(handles.gui_exp_data.fig_backbone)
                clf(handles.gui_exp_data.fig_backbone)
                close(handles.gui_exp_data.fig_backbone)
            end
         end
        
        close(handles.gui_exp_data.gui_exp)
    end
    
    sys.T = handles.sys.V;
    sys.At = handles.sys.At;
    sys.A = handles.sys.A;
    sys.f = handles.sys.f;
    sys.spv = handles.sys.spv;
    sys.symbolic = handles.sys.symbolic;
    sys.modal = handles.sys.modal;
    sys.complex = handles.sys.complex;
    sys.conservative = handles.sys.conservative; 
    sys.complex_cor = handles.complex_cor;
    sys.lambda = handles.sys.lambda;
    sys.lambda_select = handles.sys.lambda_s;
    sys.lambda_remain = handles.sys.lambda_p;
    n = handles.sys.n;
  
    
    handles.sys.ho_int_res = 0;
    

        sys.Tmodal = handles.sys.Vmodal;
        sys.sigma = handles.sys.sigma;
        sys.int_res = handles.sys.int_res;
        
        if ~handles.sys.conservative
            if n <= sys.sigma
                mode = struct('WindowStyle','nonmodal','Interpreter','tex');
                txt = 'Chosen order of SSM expansion is lower than the spectral quotient plus one.';
                h = errordlg(txt,'Choose higer order of expansion',mode);
                Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
               return 
            end
        end
        
        if handles.check_higher_int_res.Value
            try   
                [handles,sys] = check_higher_res(handles,sys);
            catch ME
                rethrow(ME)
            end
        else
            if handles.sys.int_res
                sys.int_res_vec = handles.sys.int_res_vec_init;
            end
        end
    try    
        
        pp1 = gcp('nocreate');
        if isempty(pp1)
            h = helpdlg('Starting parallel pool for the first time and detecting number of available cores.', 'Info');
            disp('----------------------------------------------------------------')
            disp('Starting parallel pool for the first time and detecting number')
            disp('of available cores.')
            disp('----------------------------------------------------------------')
            defaultProfile = parallel.defaultClusterProfile;
            myCluster = parcluster(defaultProfile);
            parpool(myCluster);
            pp2 = gcp('nocreate');
            cpuNum =pp2.NumWorkers;
            save('cluster_info.mat','cpuNum')
            
            if isvalid(h)
                close(h)
            end
        end
        

    catch ME
        rethrow(ME)
    end
        
    try
        [SSM_proj] = compute_SSM(sys,n);
    catch ME
        Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
        rethrow(ME)
    end
    
    if  handles.sys.complex
        set(handles.push_update_z,'Enable','Off')
        set(handles.push_auto,'Enable','Off')
        set(handles.push_update,'Enable','Off')
        set(handles.push_update_domain,'Enable','Off')
        set(handles.push_openint,'Enable','Off')
        set(handles.open_SSM_expr,'Enable','On')
        set(handles.open_invar,'Enable','Off')
    else
        
        try
            handles = plot_SSM(handles);
        catch MA
            Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
            rethrow(MA)
        end
        
        set(handles.push_update_z,'Enable','On')
        set(handles.push_auto,'Enable','On')
        set(handles.push_update,'Enable','On')
        set(handles.push_update_domain,'Enable','On')
        set(handles.push_openint,'Enable','On')
        set(handles.open_SSM_expr,'Enable','On')
        
        if ~handles.sys.conservative
            set(handles.open_invar,'Enable','On')
        else
            set(handles.open_invar,'Enable','Off')
        end
    end
    
    Spinner(handles.panel_spinner,'Done','Stop',handles.spinner.position);
    guidata(hObject, handles)

function input_nonlinear_order_Callback(hObject, eventdata, handles)
% hObject    handle to input_nonlinear_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_nonlinear_order as text
%        str2double(get(hObject,'String')) returns contents of input_nonlinear_order as a double


% --- Executes during object creation, after setting all properties.
function input_nonlinear_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_nonlinear_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_order.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popup_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_order contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_order


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_update.
function push_update_Callback(hObject, eventdata, handles)
% hObject    handle to push_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.plot.xmin = str2num(get(handles.plot_xmin,'String'));
    handles.plot.xmax = str2num(get(handles.plot_xmax,'String')); 
    handles.plot.ymin = str2num(get(handles.plot_ymin,'String'));
    handles.plot.ymax = str2num(get(handles.plot_ymax,'String'));
    handles.plot.zmin = str2num(get(handles.plot_zmin,'String'));
    handles.plot.zmax = str2num(get(handles.plot_zmax,'String'));
    
    
    if isempty(handles.plot.xmin) || isempty(handles.plot.xmax) || isempty(handles.plot.ymin) || isempty(handles.plot.ymax) || isempty(handles.plot.zmin) || isempty(handles.plot.zmax)
        mode = struct('WindowStyle','nonmodal','Interpreter','tex');
        txt = 'Please enter an upper and lower bound for the plot range';
        h = errordlg(txt,'Input error',mode);
        return;
    end

    set(handles.axis_handle,'XLim',[handles.plot.xmin handles.plot.xmax],...
                            'YLim',[handles.plot.ymin handles.plot.ymax],...
                            'ZLim',[handles.plot.zmin handles.plot.zmax]);
          
    % Update handles structure
    guidata(hObject, handles);





function plot_r_max_Callback(hObject, eventdata, handles)
% hObject    handle to plot_r_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_r_max as text
%        str2double(get(hObject,'String')) returns contents of plot_r_max as a double


% --- Executes during object creation, after setting all properties.
function plot_r_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_r_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_z1min_Callback(hObject, eventdata, handles)
% hObject    handle to plot_z1min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_z1min as text
%        str2double(get(hObject,'String')) returns contents of plot_z1min as a double


% --- Executes during object creation, after setting all properties.
function plot_z1min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_z1min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_z2max_Callback(hObject, eventdata, handles)
% hObject    handle to plot_z2max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_z2max as text
%        str2double(get(hObject,'String')) returns contents of plot_z2max as a double


% --- Executes during object creation, after setting all properties.
function plot_z2max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_z2max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_z2min_Callback(hObject, eventdata, handles)
% hObject    handle to plot_z2min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_z2min as text
%        str2double(get(hObject,'String')) returns contents of plot_z2min as a double


% --- Executes during object creation, after setting all properties.
function plot_z2min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_z2min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_xmax_Callback(hObject, eventdata, handles)
% hObject    handle to plot_xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_xmax as text
%        str2double(get(hObject,'String')) returns contents of plot_xmax as a double


% --- Executes during object creation, after setting all properties.
function plot_xmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_xmin_Callback(hObject, eventdata, handles)
% hObject    handle to plot_xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_xmin as text
%        str2double(get(hObject,'String')) returns contents of plot_xmin as a double


% --- Executes during object creation, after setting all properties.
function plot_xmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_ymax_Callback(hObject, eventdata, handles)
% hObject    handle to plot_ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_ymax as text
%        str2double(get(hObject,'String')) returns contents of plot_ymax as a double


% --- Executes during object creation, after setting all properties.
function plot_ymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_ymin_Callback(hObject, eventdata, handles)
% hObject    handle to plot_ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_ymin as text
%        str2double(get(hObject,'String')) returns contents of plot_ymin as a double


% --- Executes during object creation, after setting all properties.
function plot_ymin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_zmax_Callback(hObject, eventdata, handles)
% hObject    handle to plot_zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_zmax as text
%        str2double(get(hObject,'String')) returns contents of plot_zmax as a double


% --- Executes during object creation, after setting all properties.
function plot_zmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_zmin_Callback(hObject, eventdata, handles)
% hObject    handle to plot_zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_zmin as text
%        str2double(get(hObject,'String')) returns contents of plot_zmin as a double


% --- Executes during object creation, after setting all properties.
function plot_zmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_auto.
function push_auto_Callback(hObject, eventdata, handles)
% hObject    handle to push_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.axis_handle,'XlimMode','auto')
    set(handles.axis_handle,'YlimMode','auto')
    set(handles.axis_handle,'ZlimMode','auto')


% --- Executes on button press in push_update_domain.
function push_update_domain_Callback(hObject, eventdata, handles)
% hObject    handle to push_update_domain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    Spinner(handles.panel_spinner,'Busy','Start',handles.spinner.position);
    handles.plot.r_max = str2num(get(handles.plot_r_max,'String'));
    
    if isempty(handles.plot.r_max) 
        mode = struct('WindowStyle','nonmodal','Interpreter','tex');
        txt = 'Please enter a numerical value for the parameterization domain';
        h = errordlg(txt,'Input error',mode);
        Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
        return
    end
    
    handles = plot_SSM(handles);
    Spinner(handles.panel_spinner,'Done','Stop',handles.spinner.position);
    guidata(hObject, handles)


% --- Executes on selection change in popup_z_cor.
function popup_z_cor_Callback(hObject, eventdata, handles)
% hObject    handle to popup_z_cor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(handles.popup_z_cor,'Value');  
handles.plot.select(3) = val; 

% Update handles structure
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popup_z_cor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_z_cor


% --- Executes during object creation, after setting all properties.
function popup_z_cor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_z_cor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_update_z.
function push_update_z_Callback(hObject, eventdata, handles)
% hObject    handle to push_update_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = plot_SSM(handles);
guidata(hObject, handles)
% plot_SSM(hObject, eventdata, handles)


% --- Executes on button press in push_openint.
function push_openint_Callback(hObject, eventdata, handles)
% hObject    handle to push_openint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Spinner(handles.panel_spinner,'Busy','Start',handles.spinner.position);
SSM_int
Spinner(handles.panel_spinner,'Done','Stop',handles.spinner.position);


% --- Executes on button press in push_loadmech.
function push_loadmech_Callback(hObject, eventdata, handles)
% hObject    handle to push_loadmech (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Spinner(handles.panel_spinner,'Busy','Start',handles.spinner.position);
if isfield(handles,'load_ex')
    
    set(handles.push_loadmech,'Enable','Off')
    
    if handles.load_ex == 1
        struct = load(fullfile(pwd,'Examples','2DOF_inner_res.mat'));
        FileName = '2DOF Near-Inner Resonance';
    elseif handles.load_ex == 2
        struct = load(fullfile(pwd,'Examples','2DOF_outer_res.mat'));
        FileName = '2DOF Near-Outer Resonance';
    elseif handles.load_ex == 3
        struct = load(fullfile(pwd,'Examples','Beam.mat'));
        FileName = 'Nonlinear Damped Timoshenko Beam';
    end
    
    handles = rmfield(handles,'load_ex');
    set(handles.push_loadmech,'Enable','On')
    guidata(hObject, handles);
    
else
    [FileName,PathName] = uigetfile('*.mat','Select mat file');
    if FileName==0, return, end
    struct = load(fullfile(PathName,FileName));  
    set(handles.popup_load_example,'Value',1);
end

if isfield(struct,'M') && isfield(struct,'C')  && isfield(struct,'K') && isfield(struct,'f')  && isfield(struct,'x') && isfield(struct,'xd') && isfield(struct,'scaling') && isfield(struct,'conservative')
    handles.sys.M = struct.M;
    handles.sys.C = struct.C;
    handles.sys.K = struct.K;
    spv = [struct.x;struct.xd];
    handles.sys.spv = spv;
    handles.sys.scaling = struct.scaling;
    handles.sys.conservative = struct.conservative;
    
    xstring = sym2char(struct.x,[numel(struct.x),1]);
    xdstring = sym2char(struct.xd,[numel(struct.x),1]);
    
    Mstring = sym2char(struct.M,[numel(struct.x),numel(struct.x)]);
    Cstring = sym2char(struct.C,[numel(struct.x),numel(struct.x)]);
    Kstring = sym2char(struct.K,[numel(struct.x),numel(struct.x)]);
    fstring = sym2char(struct.f,[numel(struct.x),1]);
    scalingstring = sprintf('%d',struct.scaling);
        
    set(handles.input_pos,'String',xstring);
    set(handles.input_vel,'String',xdstring);
    set(handles.input_mass,'String',Mstring);
    set(handles.input_damping,'String',Cstring);
    set(handles.input_stiffness,'String',Kstring);
    set(handles.input_nonlinear,'String',fstring);
    set(handles.scaling_factor,'String',scalingstring);
    
    if handles.sys.conservative
         set(handles.radio_con_yes,'Value',1);
         set(handles.input_damping,'Enable','Off');
    else
         set(handles.radio_con_no,'Value',1);
         set(handles.input_damping,'Enable','On');
    end 
    
    set(handles.radio_modal,'Enable','On');

    clear_spectral_handles(hObject, eventdata, handles)
    set(handles.text_file_loaded,'String', sprintf('MAT-file loaded: %s',FileName));
    set(handles.text_file_loaded,'Visible','On');  
    Spinner(handles.panel_spinner,'Done','Stop',handles.spinner.position);
    % Update handles structure
    guidata(hObject, handles);
    
else
    errordlg('Error loading input file','Loading Error');
    Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
end



function clear_spectral_handles(hObject, eventdata, handles)
    
    %-- Reinitialize GUI --%
    set(handles.push_clear_loaded,'Enable','Off');
    set(handles.text_file_loaded,'String', '')
    set(handles.text_file_loaded,'Visible','Off')
    set(handles.list_lambda,'String','');
    set(handles.text_nonres,'String','');
    set(handles.text_nonres_2,'String','');
    set(handles.text_sigma,'String','');
    set(handles.text_order,'String','5');
    set(handles.text_order,'Enable','Off');
    set(handles.push_select_lambda,'Enable','Off');
   
    handles.plot.select = 1:3;
    set(handles.popup_x_cor,'Value',1);
    set(handles.popup_y_cor,'Value',1);
    set(handles.popup_z_cor,'Value',1);
    set(handles.popup_x_cor,'String',{'Choose'});
    set(handles.popup_y_cor,'String',{'Choose'});
    set(handles.popup_z_cor,'String',{'Choose'});
    
    handles.compute.lambda = 0;
    handles.compute.order = 1;
    check_status_compute(hObject, eventdata, handles)
    
    set(handles.push_update_z,'Enable','Off')
    set(handles.push_auto,'Enable','Off')
    set(handles.push_update,'Enable','Off')
    set(handles.push_update_domain,'Enable','Off')
    set(handles.push_openint,'Enable','Off')
    set(handles.open_SSM_expr,'Enable','Off')
    set(handles.open_invar,'Enable','Off')
    set(handles.check_higher_int_res,'Enable','Off');
    set(handles.check_higher_int_res,'Value',1);
    
    %-- Close figures and open GUI windows --%
    if ~isempty(handles.fig)
        if isvalid(handles.fig)
            clf(handles.fig)
            close(handles.fig)
        end
    end
    
    h1 = findobj('Tag','gui_plot');
    h2 = findobj('Tag','gui_exp');
    h3 = findobj('Tag','gui_invar');
    
    if ~isempty(h1)
        handles.gui_plot_data = guidata(h1);
        close(handles.gui_plot_data.gui_plot)
    end
    
    if ~isempty(h2)
        handles.gui_exp_data = guidata(h2);
         if ~isempty(handles.gui_exp_data.fig_backbone)
            if isvalid(handles.gui_exp_data.fig_backbone)
                clf(handles.gui_exp_data.fig_backbone)
                close(handles.gui_exp_data.fig_backbone)
            end
         end      
        close(handles.gui_exp_data.gui_exp)
    end
    
    if ~isempty(h3)
        handles.gui_invar_data = guidata(h3);
        close(handles.gui_invar_data.gui_invar)
    end
    Spinner(handles.panel_spinner,'Done','Stop',handles.spinner.position);
    % Update handles structure
    guidata(hObject, handles);
    
    
    


% --- Executes on button press in push_clear_loaded.
function push_clear_loaded_Callback(hObject, eventdata, handles)
% hObject    handle to push_clear_loaded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
    
    clear_spectral_handles(hObject, eventdata, handles)
    handles = rmfield(handles,'sys');

    % Update handles structure
    guidata(hObject, handles);
    


% --- Executes on button press in radio_no.
function radio_no_Callback(hObject, eventdata, handles)
% hObject    handle to radio_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_no

    if get(hObject,'Value') 
        set(handles.input_parv,'Enable','Off')
    end

    set(handles.radio_modal,'Enable','On');
    
    % Update handles structure
    guidata(hObject, handles);




function input_parv_Callback(hObject, eventdata, handles)
% hObject    handle to input_parv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_parv as text
%        str2double(get(hObject,'String')) returns contents of input_parv as a double


% --- Executes during object creation, after setting all properties.
function input_parv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_parv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_yes.
function radio_yes_Callback(hObject, eventdata, handles)
% hObject    handle to radio_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_yes

    if get(hObject,'Value') 
        set(handles.input_parv,'Enable','On')
    end
    
    set(handles.radio_modal,'Enable','Off')
    set(handles.radio_physical,'Value',1)
    
    set(handles.text_nonres,'String','');
    set(handles.text_nonres_2,'String','');
    set(handles.text_sigma,'String','');

    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in open_SSM_expr.
function open_SSM_expr_Callback(hObject, eventdata, handles)
% hObject    handle to open_SSM_expr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Spinner(handles.panel_spinner,'Busy','Start',handles.spinner.position);
SSM_exp
Spinner(handles.panel_spinner,'Done','Stop',handles.spinner.position);


% --- Executes on selection change in popup_x_cor.
function popup_x_cor_Callback(hObject, eventdata, handles)
% hObject    handle to popup_x_cor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_x_cor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_x_cor

    val = get(handles.popup_x_cor,'Value');
    handles.plot.select(1) = val; 

    % Update handles structure
    guidata(hObject, handles);

    


% --- Executes during object creation, after setting all properties.
function popup_x_cor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_x_cor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_y_cor.
function popup_y_cor_Callback(hObject, eventdata, handles)
% hObject    handle to popup_y_cor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_y_cor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_y_cor

    val = get(handles.popup_y_cor,'Value');
    handles.plot.select(2) = val; 

    % Update handles structure
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_y_cor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_y_cor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_physical.
function radio_physical_Callback(hObject, eventdata, handles)
% hObject    handle to radio_physical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_physical


% --- Executes on button press in radio_modal.
function radio_modal_Callback(hObject, eventdata, handles)
% hObject    handle to radio_modal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_modal


% --- Executes on button press in push_save.
function push_save_Callback(hObject, eventdata, handles)
% hObject    handle to push_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

            
 radio_select_con = get(handles.buttongroup_con,'SelectedObject');
    
    switch radio_select_con.String
        case 'Yes'
             handles.sys.conservative = 1;
        case 'No'
             handles.sys.conservative = 0;
    end  
    
    
    M_string = get(handles.input_mass,'String'); 
    C_string = get(handles.input_damping,'String');
    K_string = get(handles.input_stiffness,'String'); 

    q_state = get(handles.input_pos,'String');
    qd_state = get(handles.input_vel,'String');

    q_state(q_state==' ') = '';
    q_state([strfind(q_state,'['),strfind(q_state,']')]) = [];
    q = sym(strread(q_state,'%s','delimiter',';'));

    qd_state(qd_state==' ') = '';
    qd_state([strfind(qd_state,'['),strfind(qd_state,']')]) = [];
    qd = sym(strread(qd_state,'%s','delimiter',';'));
    
    scaling = str2double(get(handles.scaling_factor,'String'));
    

    if isempty(q_state) || isempty(qd_state)
        errordlg('State variables are not specified','Input Error');
    elseif isempty(M_string) || isempty(C_string) || isempty(K_string)
        errordlg('System matrices are not specified','Input Error');
    elseif isempty(get(handles.input_nonlinear,'String'))
        errordlg('Nonlinear vector is not specified','Input Error');
    elseif isempty(get(handles.scaling_factor,'String'))
        errordlg('Scaling vector is not specified','Input Error');
    else
               
    spv  = [q;qd];
    syms(spv(:))
     
    f = eval(get(handles.input_nonlinear,'String'));
        
    M =  sym(M_string);
    C =  sym(C_string);
    K =  sym(K_string);
    x = q;
    xd = qd;
    conservative = handles.sys.conservative;

    [FileName,PathName] = uiputfile('*.mat','Save mat file');
    if FileName==0, return, end
        save(sprintf('%s%s',PathName,FileName),'M','C','K','f','x','xd','scaling','conservative');
    end
 

% --- Executes on mouse press over figure background.
function gui_main_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to gui_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radio_con_yes.
function radio_con_yes_Callback(hObject, eventdata, handles)
% hObject    handle to radio_con_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_con_yes

    if get(hObject,'Value') 
        set(handles.input_damping,'Enable','Off')
    end
    
    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in radio_con_no.
function radio_con_no_Callback(hObject, eventdata, handles)
% hObject    handle to radio_con_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_con_no

  if get(hObject,'Value') 
        set(handles.input_damping,'Enable','On')
  end


% --- Executes when panel_url is resized.
function panel_url_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to panel_url (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when gui_main is resized.
function gui_main_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to gui_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in pop_academic_ex.
function pop_academic_ex_Callback(hObject, eventdata, handles)
% hObject    handle to pop_academic_ex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_academic_ex contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_academic_ex


% --- Executes during object creation, after setting all properties.
function pop_academic_ex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_academic_ex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_load_example.
function popup_load_example_Callback(hObject, eventdata, handles)
% hObject    handle to popup_load_example (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


str = get(hObject, 'String');
val = get(hObject,'Value');

switch str{val}
    
case '2DOF Near-Inner Resonance' 
    
    handles.load_ex = 1;
    push_loadmech_Callback(hObject, eventdata, handles)
    handles = guidata(hObject); 
case '2DOF Near-Outer Resonance'
    
    handles.load_ex = 2;
    push_loadmech_Callback(hObject, eventdata, handles)
    handles = guidata(hObject); 
case 'Nonlinear Damped Timoshenko Beam' 
    
    handles.load_ex = 3;
    push_loadmech_Callback(hObject, eventdata, handles)
    handles = guidata(hObject); 
otherwise
end

% Save the handles structure.
guidata(hObject,handles)



% Hints: contents = cellstr(get(hObject,'String')) returns popup_load_example contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_load_example


% --- Executes during object creation, after setting all properties.
function popup_load_example_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_load_example (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close gui_main.
function gui_main_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to gui_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

     if ~isempty(handles.fig)
        if isvalid(handles.fig)
            clf(handles.fig)
            close(handles.fig)
        end
     end
    
    h1 = findobj('Tag','gui_plot');
    h2 = findobj('Tag','gui_exp');
    
    if ~isempty(h1)
        handles.gui_plot_data = guidata(h1);
        close(handles.gui_plot_data.gui_plot)
    end
    
    if ~isempty(h2)
        handles.gui_exp_data = guidata(h2);
        
         if ~isempty(handles.gui_exp_data.fig_backbone)
            if isvalid(handles.gui_exp_data.fig_backbone)
                clf(handles.gui_exp_data.fig_backbone)
                close(handles.gui_exp_data.fig_backbone)
            end
         end

        close(handles.gui_exp_data.gui_exp)
    end
    clc
    
    
    disp('----------------------------------------------------------------')
    disp('                       GUI Closed                               ')
    disp('----------------------------------------------------------------')
 
% Hint: delete(hObject) closes the figure
delete(hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over label_logo.
function label_logo_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to label_logo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text_url.
function text_url_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text_url (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
url = 'http://www.georgehaller.com';
[stat,h] = web(url);


% --- Executes on button press in check_higher_int_res.
function check_higher_int_res_Callback(hObject, eventdata, handles)
% hObject    handle to check_higher_int_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if ~get(hObject,'Value')
%     handles.sys.int_res_vec = handles.sys.int_res_vec_init;
% end
% 
% guidata(hObject,handles)

% Hint: get(hObject,'Value') returns toggle state of check_higher_int_res



function scaling_factor_Callback(hObject, eventdata, handles)
% hObject    handle to scaling_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scaling_factor as text
%        str2double(get(hObject,'String')) returns contents of scaling_factor as a double


% --- Executes during object creation, after setting all properties.
function scaling_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scaling_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open_invar.
function open_invar_Callback(hObject, eventdata, handles)
% hObject    handle to open_invar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Spinner(handles.panel_spinner,'Busy','Start',handles.spinner.position);
SSM_invar
Spinner(handles.panel_spinner,'Done','Stop',handles.spinner.position);



function text_order_Callback(hObject, eventdata, handles)
% hObject    handle to text_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_order as text
%        str2double(get(hObject,'String')) returns contents of text_order as a double

check_status_compute(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function text_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on text_order and none of its controls.
function text_order_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to text_order (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
 
