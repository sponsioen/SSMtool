function [output] = compute_SSM(sys,n)
clc
disp('----------------------------------------------------------------')
disp('Automated computation of Spectral Submanifolds for autonomous')
disp('mechanical systems V1.0 (2017)')
disp(' ')
disp('Sten Ponsioen (stenp@ethz.ch)')
disp('George Haller (georgehaller@ethz.ch)')
disp('----------------------------------------------------------------')

order = n;
lambda_order = sys.lambda.num([sys.lambda_select,sys.lambda_remain]);
T = sys.T;
At = sys.At;
f = sym(zeros(numel(sys.spv),1));
check_eq = abs(sum(subs(sys.f,sys.spv,zeros(numel(sys.spv),1))));
mode = struct('WindowStyle','nonmodal','Interpreter','tex');

if check_eq > 1e-10
    errordlg('The equilibrium point of the system is non-trivial','Non-trivial equilibrium point', mode); 
    return;
end
  
for i = 1:numel(sys.spv)
   f(i) = simplify(taylor(sys.f(i),sys.spv.','ExpansionPoint',0,'Order',2)-taylor(sys.f(i),sys.spv.','ExpansionPoint',0,'Order',1));
end

if numel(find(f==0))<numel(sys.spv)
    errordlg('The Taylor expansion of the nonlinear force vector contains linear terms.','Order Error', mode); 
    return;
end

for i = 1:numel(sys.spv)
    f(i) = simplify(taylor(sys.f(i),sys.spv.','ExpansionPoint',0,'Order',order+1));
end

f = expand(f);
spv = sys.spv;
modal = sys.modal;
complex = sys.complex;
int_res = sys.int_res;
Tmodal = sys.Tmodal;
sys_dim = numel(spv);
syms z1 z2
X1 = [z1; z2];
R1 = diag(sys.lambda.num(sys.lambda_select));
K1 = [eye(2);zeros(sys_dim-2,2)];
K = cell(1,order); R = cell(1,order); X = cell(1,order);
K{1} = K1; R{1} = R1; X{1} = X1;

for i = 2:order
    K{i} = zeros(sys_dim,2^i);
    R{i} = zeros(2,2^i);
    c = combinator(2,i,'p','r');
    b1 = sum(c==1,2);
    b2 = sum(c==2,2);
    X{i} = z1.^b1.*z2.^b2;
end

try
    G = matGV2(f,spv,T,order);
catch GR
    Spinner(handles.panel_spinner,'Error','Stop',handles.spinner.position);
    rethrow(GR)
end
RI = {eye(2),R1};
hwait=waitbar(0,'Computing Spectral Submanifold'); 

disp(' ')
disp('----------------   Locating internal resonances   ---------------- ')
disp(' ')
    if int_res
        disp('Resonant terms:')
                loc_R = cell(2,order);
                int_res_orders = unique(sum(cat(1,sys.int_res_vec{:,:}),2));
                p_r = cell(1,max(int_res_orders));
                disp('z1^i*z2^j:')
                disp(cat(1,sys.int_res_vec{:,:}))
                for j = int_res_orders.'
                    p_r{j} = combinator(numel(X1),j,'p','r');
                end
                for q=1:2 
                    coef_int_q = cat(1,sys.int_res_vec{q,:});
                    coef_int_q_orders = sum(coef_int_q,2);
                    for k=1:size(coef_int_q,1)
                        int_vec = sort([repmat(1,[1,coef_int_q(k,1)]),repmat(2,[1,coef_int_q(k,2)])],2);
                        loc  = find(any(all(bsxfun(@eq,int_vec,sort(p_r{coef_int_q_orders(k)},2)),2),3)).';
                        loc_R{q,coef_int_q_orders(k)} = [loc_R{q,coef_int_q_orders(k)},loc];
                    end
                end
        disp(' ')
        disp('Done.')
    else
       disp('No resonant terms found.')
    end

disp(' ')
disp('------------------   Computing W(z) and R(z)   ------------------- ')
disp(' ')
for n=2:order
    tic
    fprintf('Solving order: %d \n',n)
    hwait=updatewaitbar(hwait,n,order);
    free_dofs = 1:sys_dim*2^n;
    c = combinator(2,n,'p','r');
    RnI = sum(lambda_order(c),2);
    At_vec = repmat(diag(At),[1,2^n]);
    Ae_vec = repmat(RnI.',[sys_dim,1]);
    LHS_n = reshape(At_vec-Ae_vec,[],1);  
    KR = zeros(sys_dim,2^n);
    G_check = G{n};
    G_check_bool = numel(G_check(:))>1;
    if  G_check_bool
        GK1 =  kronGK1n(n,K,G);
    end
    GK = zeros(sys_dim,2^n);
    for m = 2:n-1
        KR = KR + kronKR(n,m,K,R);
        GK = GK + kronGK(n,m,K,G);  
    end

    if int_res
        Rn_dummy = R{n};
        KR_dummy = KR(1:2,:);
        if  G_check_bool
            GK1_dummy = GK1(1:2,:);      
        end
        GK_dummy = GK(1:2,:);
        for q = 1:2
            if numel(loc_R{q,n}) > 0                  
                if G_check_bool 
                    Rn_dummy(q,loc_R{q,n}) = -KR_dummy(q,loc_R{q,n}) + GK1_dummy(q,loc_R{q,n}) + GK_dummy(q,loc_R{q,n});
                else 
                    Rn_dummy(q,loc_R{q,n}) = -KR_dummy(q,loc_R{q,n}) + GK_dummy(q,loc_R{q,n});
                end     
            end
        end
        R{n} = Rn_dummy;
        dof = 1:sys_dim*2^n;
        constraints_u =  (loc_R{1,n}-1)*sys_dim+1;
        constraints_l =  (loc_R{2,n}-1)*sys_dim+2;
        constraints = union(constraints_u, constraints_l);
        free_dofs = setdiff(dof,constraints);
    end

    if G_check_bool
        RHS_n = K1*R{n} + KR - GK1 - GK;  
    else
        RHS_n = K1*R{n} + KR - GK;   
    end
    RHS_n_vec = reshape(RHS_n,[sys_dim*2^n,1]);
    K_n_vec_full = zeros(sys_dim*2^n,1);
    K_diag = LHS_n(free_dofs);
    P_vec = RHS_n_vec(free_dofs);
    
    Np = numel(K_diag);
    cluster = load('cluster_info.mat');
    cpuNum = cluster.cpuNum;
    id = round(linspace(0,Np,cpuNum+1));
    spmd 
        range = id(labindex)+ 1:id(labindex+1);  
        sol_con = P_vec(range)./K_diag(range);   
    end
    K_n_vec_full(free_dofs) = cat(1,sol_con{:});
    K_n = reshape(K_n_vec_full,[sys_dim,2^n]);
    K{n} = K_n;
    toc
    disp(' ')
end
disp('Done.')
close(hwait)

SSM_FPS = zeros(sys_dim,1);
R_FPS = zeros(2,1);

for i = 1:order
    SSM_FPS = SSM_FPS + round(K{i},17)*X{i};
    R_FPS = R_FPS + round(R{i},17)*X{i};
end

if modal
    SSM_proj = Tmodal\T*SSM_FPS;
elseif complex
    SSM_proj = SSM_FPS;
else
    SSM_proj = T*SSM_FPS;
end

input_vars_SSM = cell(1,numel(SSM_proj));
output_vars_SSM = cell(1,numel(SSM_proj));
input_vars_R = cell(1,numel(R_FPS));
output_vars_R = cell(1,numel(R_FPS));

for i = 1:numel(SSM_proj)
input_vars_SSM{i} = SSM_proj(i);
output_vars_SSM{i} = strcat('q_', num2str(i));
end

for i = 1:numel(R_FPS)
input_vars_R{i} = R_FPS(i);
output_vars_R{i} = strcat('zd_', num2str(i));
end

disp(' ')
disp('------   Writing SSM and reduced dynamics to function file   ----- ')
disp(' ')

nm = load('cs.mat');
mkdir('./Data/',strcat('run_',nm.folder_id))
cd(strcat('./Data/run_',nm.folder_id))
matlabFunction(input_vars_SSM{:},'file',strcat('SSM_function_',nm.folder_id),'Vars',X1,'Outputs',output_vars_SSM);
save(strcat('sys_',nm.folder_id),'sys')
disp('W(z) written to: SSM_function.m')
matlabFunction(input_vars_R{:},'file',strcat('R_function_',nm.folder_id),'Vars',X1,'Outputs',output_vars_R);
disp('R(z) written to: R_function.m')
cd('../..')
addpath(genpath(strcat(pwd,'/Data')));
matlabFunction(input_vars_R{:},'file','R_sub_function','Vars',X1,'Outputs',output_vars_R);

disp(' ')
disp('Done.')
output = SSM_proj;

end

