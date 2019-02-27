function [te,ystate,Rm,err,tevent_out] = int_red_dyn_em(t_end,eps_init,y0,set_options,nm,sys)

SSM_function = @(z1,z2)eval(strcat('SSM_function_',nm.folder_id,'(z1,z2)'));
maxstep = set_options.maxstep;
complex_cor = set_options.complex_cor;
options = odeset('RelTol',1e-11,'AbsTol',1e-11,'MaxStep',maxstep,'Events',@event_function);
err = zeros(1,size(y0,1));
tevent_out = cell(1,size(y0,1));
te = cell(1,size(y0,1));
ystate = cell(1,size(y0,1));
Rm = cell(1,size(y0,1));

for i=1:size(y0,1)
    ythres = eps_init;
    [tInter,ystateInter,tevent,yevent,ievent] = ode45(@(t,y) odefun(t,y),[0,t_end],y0(i,:),options);

    if isempty(ievent)
        err(i) = 1;
        return;
    end

    tlin = linspace(0,tevent,set_options.interp);
    ystateInter = interp1(tInter,ystateInter,tlin,'spline');
    te{i} = tlin;
    ystate{i} = ystateInter;
    tevent_out{i} = tevent;
    R = cell(1,numel(sys.spv));
    Rho = ystateInter(:,1);
    Theta = ystateInter(:,2);

    if complex_cor
        [R{:}] = SSM_function(Rho.*exp(1i.*Theta),Rho.*exp(-1i.*Theta));
    else
        [R{:}] = SSM_function(Rho.*cos(Theta),Rho.*sin(Theta));
    end    

    Rm{i} = real([R{:}]);

end

function dq = odefun(t,y)
    if complex_cor
        [zd1,zd2] = R_sub_EM_function(y(1)*exp(1i*y(2)),y(1)*exp(-1i.*y(2)));
        alpha = real(zd1);
        beta = imag(zd1);
    else
        [zd1,zd2] = R_sub_EM_function(y(1)*cos(y(2)),y(1)*sin(y(2)));
        alpha = real(zd1);
        beta = real(zd2);
    end

    qd1 = alpha*cos(y(2))+ beta*sin(y(2));
    qd2 = (1/y(1))*(beta*cos(y(2))-alpha*sin(y(2)));
    dq = [qd1;qd2];
end

function [value,isterminal,direction] = event_function(t,y)
    value = y(1) - ythres; 
    isterminal = 1; 
    direction = -1; 
end

end