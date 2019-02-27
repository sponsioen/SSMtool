function [te,ystate] = int_red_dyn(t_end,y0,Z,set_options)

maxstep = set_options.maxstep;
complex_cor = set_options.complex_cor;
nm = load('cs.mat');
R = cell(2,1);
[R{:}] = eval(strcat('R_function_',nm.folder_id,'(Z(1),Z(2))'));
matlabFunction(R{:},'file','R_sub_function','Vars',Z);
options = odeset('RelTol',1e-11,'AbsTol',1e-15,'MaxStep',maxstep);
te = cell(1,size(y0,1));
ystate = cell(1,size(y0,1));

for i=1:size(y0,1)
[tInter,ystateInter] = ode45(@(t,y) odefun(t,y),[0,t_end],y0(i,:),options);
te{i} = tInter(1:end,:);
ystate{i} = ystateInter(1:end,:);
end

te = cat(1,te{:});
ystate = cat(1,ystate{:});

function dq = odefun(t,y)
    if complex_cor
        [zd1,zd2] = R_sub_function(y(1)*exp(1i*y(2)),y(1)*exp(-1i.*y(2)));
        alpha = real(zd1);
        beta = imag(zd1); 
    else
        [zd1,zd2] = R_sub_function(y(1)*cos(y(2)),y(1)*sin(y(2)));
        alpha = real(zd1);
        beta = real(zd2);
    end
    qd1 = alpha*cos(y(2))+ beta*sin(y(2));
    qd2 = (1/y(1))*(beta*cos(y(2))-alpha*sin(y(2)));
    dq = [qd1;qd2];
end

end