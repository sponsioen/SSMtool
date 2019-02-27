function [te,ystate] = int_dyn(sys,t_end,y0,set_options)

maxstep = set_options.maxstep;
options = odeset('RelTol',1e-11,'AbsTol',1e-11,'MaxStep',maxstep);

A = sys.A;
f = sys.f;
spv = sys.spv;

dyn_sys = A*spv + f;
matlabFunction(dyn_sys,'file','system_function','Vars',[spv]);

te = cell(1,size(y0,1));
ystate = cell(1,size(y0,1));

for i=1:size(y0,1)
[tInter,ystateInter] = ode15s(@(t,y) odefun(t,y),[0,t_end],y0(i,:),options);
te{i} = tInter(1:end,:);
ystate{i} = ystateInter(1:end,:);
end

te = cat(1,te{:});
ystate = cat(1,ystate{:});

function dq = odefun(t,y)
      state = num2cell(y.');
      dq = system_function(state{:});
end

end