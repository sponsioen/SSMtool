function [handles_out] = orderT(A,handles)

getValues = handles.sys.lambda_s; 
p = handles.sys.lambda_p;
T = zeros(size(A));
Tmodal = zeros(size(A));
lambda_select = handles.sys.lambda_select;
lambda = handles.sys.lambda.sym;

if (double(abs(imag(lambda_select(1)))) < 1e-8) && (double(abs(imag(lambda_select(2)))) < 1e-8)
    handles.complex_cor = 0;
    T(:,1)  = handles.sys.T(:,getValues(1));
    T(:,2)  = handles.sys.T(:,getValues(2));

    Tmodal(:,1)  = handles.sys.T(:,getValues(1));
    Tmodal(:,2)  = handles.sys.T(:,getValues(2));

    l = 3;
    k = 1;
   while k < numel(p)
        if (double(abs(imag(lambda(p(k))))) < 1e-8)
            T(:,l)  = handles.sys.T(:,p(k));
            Tmodal(:,l)  = handles.sys.T(:,p(k));
            l= l+1;
            k= k+1;

        elseif (double(abs(imag(lambda(p(k))))) > 1e-8) && (double(abs(lambda(p(k))-conj(lambda(p(k+1))))) < 1e-8)

            T(:,l)   = handles.sys.T(:,p(k));
            T(:,l+1) = handles.sys.T(:,p(k+1));

            Tmodal(:,l)    = real(handles.sys.T(:,p(k)));
            Tmodal(:,l+1)  = imag(handles.sys.T(:,p(k)));

            l = l+2;
            k = k+2;
        end
   end

   At = T\A*T;   
   handles.sys.At = At; 
   handles.sys.V = T;
   handles.sys.Vmodal = Tmodal; 

elseif (double(abs(imag(lambda_select(1)))) > 1e-8) && (double(abs(lambda_select(1)-conj(lambda_select(2)))) < 1e-8) %% CHECK
    handles.complex_cor = 1;
    T(:,1)  = handles.sys.T(:,getValues(1));
    T(:,2)  = handles.sys.T(:,getValues(2));

    Tmodal(:,1)  = real(handles.sys.T(:,getValues(1)));
    Tmodal(:,2)  = imag(handles.sys.T(:,getValues(1)));

    l = 3;
    k = 1;
   while k <= numel(p)
        if (double(abs(imag(lambda(p(k))))) < 1e-8)
            T(:,l)  = handles.sys.T(:,p(k));
            Tmodal(:,l)  = handles.sys.T(:,p(k));
            l= l+1;
            k= k+1;

        elseif (double(abs(imag(lambda(p(k))))) > 1e-8) && (double(abs(lambda(p(k))-conj(lambda(p(k+1))))) < 1e-8)

            T(:,l)   = handles.sys.T(:,p(k));
            T(:,l+1) = handles.sys.T(:,p(k+1));

            Tmodal(:,l)    = real(handles.sys.T(:,p(k)));
            Tmodal(:,l+1)  = imag(handles.sys.T(:,p(k)));

            l = l+2;
            k = k+2;
        end
   end

   At = T\A*T;    

   handles.sys.At = At; 
   handles.sys.V = T;
   handles.sys.Vmodal = Tmodal; 

else
    mode = struct('WindowStyle','nonmodal', 'Interpreter','tex');
    txt =  'Cannot form a 2D spectral subspace for the chosen set of eigenvalues.';
    h = errordlg(txt,'Selection Error', mode);
    handles.compute.lambda = 0;
    check_status_compute(hObject, eventdata, handles)
    error(txt)
end

handles_out = handles;
