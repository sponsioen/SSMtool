function [te,ystate] = int_dyn_em(sys,t_end,y0,set_options)

maxstep = set_options.maxstep;
options = odeset('RelTol',1e-11,'AbsTol',1e-11,'MaxStep',maxstep);
te = cell(1,size(y0,1));
ystate = cell(1,size(y0,1));

for i=1:size(y0,1)
    [tInter,ystateInter] = ode15s(@(t,y) odefun(t,y),[0,t_end(i)],y0(i,:),options);

    if sys.modal
        ystateInter = (sys.Tmodal\(ystateInter.')).';
    end

    tlin = linspace(0,t_end(i),set_options.interp);
    ystateInter = interp1(tInter,ystateInter,tlin,'spline');
    te{i} = tlin;
    ystate{i} = ystateInter;
end

function dq = odefun(t,y)
      state = num2cell(y.');
      dq = system_function(state{:});
end



end