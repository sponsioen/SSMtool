function [handles_out,sys_out] = check_higher_res(handles,sys)
getValues = handles.sys.lambda_s; 
lambda = handles.sys.lambda.sym;
n = handles.sys.n;
nonlin_vec =  cell(2,n);

syms a b lambda_1 lambda_2 lambda_k z1 z2
res = abs(([a;b;-1].'*[lambda_1;lambda_2;lambda_k])/(norm([a;b;-1])*norm([lambda_1;lambda_2;lambda_k])));
matlabFunction(res,'file','res_function','Vars',[a;b;lambda_1;lambda_2;lambda_k]);

if ~handles.sys.int_res && ~handles.sys.conservative
    ho_int_res = 0;
    for l = 2:n
        [~,ord] = nsumk(2,l) ;
        q = 0;
        for k = [getValues(1),getValues(2)]
            q = q+1;
            res_vec = double(res_function(ord(:,1),ord(:,2),lambda(getValues(1)),lambda(getValues(2)),lambda(k)));  
            in = find(res_vec < 5e-2);
            if ~isempty(in)
                pos = ord(in,:);
                for i = 1:size(pos,1)
                    nonlin_vec{q,l} = [nonlin_vec{q,l};pos(i,:)];
                end
                ho_int_res = 1;
            end
        end
    end  

    if ho_int_res 
        handles.sys.ho_int_res = 1;
        sys.int_res = 1;
        sys.int_res_vec =  nonlin_vec;
    end

elseif handles.sys.int_res && ~handles.sys.conservative
    sigma_int = handles.sys.sigma_int;
    ho_int_res = 0;

    for i = 1:size(handles.sys.int_res_vec_init,1)
        for j = 1:size(handles.sys.int_res_vec_init,2)
            nonlin_vec{i,j} =  handles.sys.int_res_vec_init{i,j};
        end
    end

    for l = sigma_int+1:n
        [~,ord] = nsumk(2,l) ;
        q = 0;
        for k = [getValues(1),getValues(2)]
            q = q+1;
            res_vec = double(res_function(ord(:,1),ord(:,2),lambda(getValues(1)),lambda(getValues(2)),lambda(k)));
            in = find(res_vec < 5e-2);
            if ~isempty(in)
                pos = ord(in,:);   
                for i = 1:size(pos,1)
                    nonlin_vec{q,l} = [nonlin_vec{q,l};pos(i,:)];
                end
                ho_int_res = 1;
            end  
        end
    end  

    if ho_int_res 
        sys.int_res_vec =  nonlin_vec;
    else
        sys.int_res_vec = handles.sys.int_res_vec_init;
    end
    
elseif  handles.sys.int_res && handles.sys.conservative   
    for l = 2:n
        [~,ord] = nsumk(2,l) ;
        q = 0;
        for k = [getValues(1),getValues(2)]
            q = q+1;
            res_vec = double(res_function(ord(:,1),ord(:,2),lambda(getValues(1)),lambda(getValues(2)),lambda(k)));
            in = find(res_vec < 5e-2);
            if ~isempty(in)
                pos = ord(in,:);
                for i = 1:size(pos,1)
                    nonlin_vec{q,l} = [nonlin_vec{q,l};pos(i,:)];
                end
            end 
        end
    end  
sys.int_res_vec =  nonlin_vec;
end 
 sys_out = sys;
 handles_out = handles;
             
