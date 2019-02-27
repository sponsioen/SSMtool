function [handles_out] = check_res(handles)

getValues = handles.sys.lambda_s; 
p = handles.sys.lambda_p;
lambda_select = handles.sys.lambda_select;
lambda_remain = handles.sys.lambda_remain;
lambda = handles.sys.lambda.sym;

syms a b lambda_1 lambda_2 lambda_k z1 z2
res = abs(([a;b;-1].'*[lambda_1;lambda_2;lambda_k])/(norm([a;b;-1])*norm([lambda_1;lambda_2;lambda_k])));
matlabFunction(res,'file','res_function','Vars',[a;b;lambda_1;lambda_2;lambda_k]);

if ~handles.sys.conservative
    sigma = double(floor(real(lambda_remain(end))/real(lambda_select(1))));
    if sigma > 50
        sigma = 5; 
        uiwait(msgbox('Large outer spectral quotient detected, sigma is set to a default value of 5.','Setting Spectral Quotient','warn'));
    end
    handles.sys.sigma = sigma;                   
    for l = 2:sigma
        [~,ord] = nsumk(2,l) ;
        for k = p
            for m = 1:size(ord,1)
                if double(res_function(ord(m,1),ord(m,2),lambda(getValues(1)),lambda(getValues(2)),lambda(k))) < 1e-4  
                    mode = struct('WindowStyle','nonmodal','Interpreter','tex');
                    txt = 'The external nonresonance conditions are violated.';
                    h = errordlg(txt,'External resonance detected',mode);
                    Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
                    error(txt)
                end
            end
        end
    end 
    sigma_int = sigma;
    handles.sys.sigma_int = sigma_int;   
    nonlin_vec = cell(2,sigma);
    int_res = 0;

    for l=2:sigma_int
        [~,ord] = nsumk(2,l) ;
        q = 0;    
        for k = [getValues(1),getValues(2)]
            q = q+1;                    
            for m = 1:size(ord,1)                             
                if double(res_function(ord(m,1),ord(m,2),lambda(getValues(1)),lambda(getValues(2)),lambda(k))) < 5e-2                                      
                    pos = ord(m,:);
                    nonlin_vec{q,l} = [nonlin_vec{q,l},pos];
                    set(handles.text_nonres_2,'String','Yes')
                    int_res = 1;
                end
            end 
        end
    end

    if ~int_res
        set(handles.text_nonres_2,'String','None')
        handles.sys.int_res = 0;
    else 
        handles.sys.int_res = 1;
        handles.sys.int_res_vec_init = nonlin_vec;
    end

    set(handles.text_nonres,'String','None')
    set(handles.text_sigma,'String',num2str(sigma))
else 
    handles.sys.sigma = NaN;
    handles.sys.int_res = 1;
    set(handles.check_higher_int_res,'Enable','Off');
    set(handles.check_higher_int_res,'Value',1);
    set(handles.text_nonres_2,'String','Yes')
    set(handles.text_nonres,'String','None')
    set(handles.text_sigma,'String','N/A')
end
             
handles_out = handles;