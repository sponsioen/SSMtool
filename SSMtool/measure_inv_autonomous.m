function [error_avg] = measure_inv_autonomous(rho_0,rho_1,N,T)
syms z1 z2
nm = load('cs.mat');
s = load(strcat('sys_',nm.folder_id));
SSM_function = @(z1,z2)eval(strcat('SSM_function_',nm.folder_id,'(z1,z2)'));
z10 = rho_0; 
eps_init = rho_1;
ntheta = N;
tend = T;
theta_lin = linspace(0,2*pi,ntheta).';
y0 = [z10.*ones(numel(theta_lin),1),theta_lin];
y0_state_per_modal = zeros(numel(theta_lin),numel(s.sys.spv));
y0_state_per = zeros(numel(theta_lin),numel(s.sys.spv));

for j =1:numel(theta_lin)
    y0_state_dum = cell(1,numel(s.sys.spv));
    if s.sys.modal
        if s.sys.complex_cor
            [y0_state_dum{:}] = SSM_function(y0(j,1)*exp(1i*y0(j,2)),y0(j,1)*exp(-1i.*y0(j,2)));
        else
            [y0_state_dum{:}] = SSM_function(y0(j,1)*cos(y0(j,2)),y0(j,1)*sin(y0(j,2)));
        end
        y0_state = real(s.sys.Tmodal*([y0_state_dum{:}].')).';
        y0_state_per_modal(j,:) = real([y0_state_dum{:}]);
        y0_state_per(j,:) = y0_state;
    elseif s.sys.complex
        if s.sys.complex_cor
            [y0_state_dum{:}] = SSM_function(y0(j,1)*exp(1i*y0(j,2)),y0(j,1)*exp(-1i.*y0(j,2)));
        else
            [y0_state_dum{:}] = SSM_function(y0(j,1)*cos(y0(j,2)),y0(j,1)*sin(y0(j,2)));
        end
        y0_state = real(s.sys.T*([y0_state_dum{:}].')).';
        y0_state_per(j,:) = y0_state;
    else
        if s.sys.complex_cor
            [y0_state_dum{:}] = SSM_function(y0(j,1)*exp(1i*y0(j,2)),y0(j,1)*exp(-1i.*y0(j,2)));
        else
            [y0_state_dum{:}] = SSM_function(y0(j,1)*cos(y0(j,2)),y0(j,1)*sin(y0(j,2)));
        end
        y0_state = real([y0_state_dum{:}]);
        y0_state_per(j,:) = y0_state;
    end
end


options.maxstep = 1e-2;
options.complex_cor = s.sys.complex_cor;
options.interp = 4000;

Zcor = [z1;z2];
R = cell(2,1);
[R{:}] = eval(strcat('R_function_',nm.folder_id,'(Zcor(1),Zcor(2))'));
matlabFunction(R{:},'file','R_sub_EM_function','Vars',Zcor);


A = s.sys.A;
f = s.sys.f;
spv = s.sys.spv;
dyn_sys = A*spv + f;
matlabFunction(dyn_sys,'file','system_function','Vars',[spv]);
Np = numel(y0_state_per(:,1));
cluster = load('cluster_info.mat');
cpuNum = cluster.cpuNum;
id = round(linspace(0,Np,cpuNum+1));

try
    spmd 
    range = id(labindex)+ 1:id(labindex+1);
    [t_red,ystate_red,Rm,err,tevent] = int_red_dyn_em(tend,eps_init,y0(range,:),options,nm,s.sys);
    end

    if sum(cat(2,err{:}))>0
        mode = struct('WindowStyle','nonmodal','Interpreter','tex');
        errordlg('Maximum integration time reached. Either increase the integration time or increase rho\_epsilon.','Error max. int. time reached', mode); 
    end
catch NE
    rethrow(NE)
end

tevent = cat(2,tevent{:});
tevent = [tevent{:}];

try
    spmd 
    range = id(labindex)+ 1:id(labindex+1);
    [t_full,ystate_full] = int_dyn_em(s.sys,tevent(range),y0_state_per(range,:),options);
    end
    
catch NE
    rethrow(NE)
end


error = zeros(options.interp,ntheta);
numlab = zeros(1,cpuNum);
for i=1:cpuNum 
   numlab(i) = numel(Rm{i});
end

k=1;
for i=1:cpuNum
 R_inter = Rm{i};
 ystate_full_inter = ystate_full{i};
    for j=1:numlab(i)
         R_inter_2 = R_inter{j};
         ystate_full_inter_2 = ystate_full_inter{j};
        
        for l =1:options.interp
            error(l,k) = sqrt((R_inter_2(l,:) - ystate_full_inter_2(l,:))*(R_inter_2(l,:)-ystate_full_inter_2(l,:)).');
        end     
        k=k+1;
    end
end

syms rho theta real
S = cell(1,numel(s.sys.spv));

if s.sys.complex_cor
    [S{:}] = SSM_function(rho.*exp(1i.*theta),rho.*exp(-1i.*theta));
else
    [S{:}] = SSM_function(rho.*cos(theta),rho.*sin(theta));
end

if s.sys.modal
    S_vec = s.sys.Tmodal*([S{:}].');
elseif s.sys.complex
    S_vec = s.sys.T*([S{:}].');
else
    S_vec = [S{:}].';
end

h = @(rho,theta)real(eval(sqrt(S_vec.'*S_vec)));
theta_linspace = linspace(0,2*pi,1000);
max_amp = max(h(z10,theta_linspace));
error_avg = mean(max(error,[],1))/max_amp;

end
