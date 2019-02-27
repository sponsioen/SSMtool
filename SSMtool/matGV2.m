function varout = matGV2(f,x,T,order)
mode = struct('WindowStyle','nonmodal','Interpreter','tex');
G = cell(1,order);
G{1} = 0;
spv = x;
hwait=waitbar(0,'Building up matrices'); 
k = 0;
c = cell(1,numel(spv));
t = cell(1,numel(spv));
t_num = zeros(numel(spv),1);
for i=1:numel(spv)
    [c{i},t{i}] = coeffs(f(i),spv.');
    t_num(i) = numel(t{i});
end
coeff_list = cat(2,c{:});
monomial_list = cat(2,t{:});
id_mon = cell(numel(monomial_list),1);
n_order = zeros(numel(monomial_list),1);
parfor i=1:numel(monomial_list)
    fact_mon = factor(monomial_list(i));
    index_mon = double(jacobian(fact_mon.',spv));
    mon_order = size(index_mon,1);
    id_mon{i} = [mon_order,sum(index_mon.*repmat(1:numel(spv),[mon_order,1]),2).']; 
    n_order(i) = numel(fact_mon);
end
nonlin_order = unique(n_order).';
iterate_remaining = 2:order;
iterate_remaining(nonlin_order-1)=[];

if max(nonlin_order) >  order 
     h = errordlg('The order of the nonlinearities is higher than the expansion order of the SSM.','Order Error', mode); 
     return;
elseif isempty(nonlin_order) 
     h = errordlg('The Taylor expansion of the nonlinear force vector for the current SSM order is zero.','Order Error', mode); 
     return;
end

disp(' ')
disp('--------   Scanning structure of nonlinear force vector  --------- ')
disp(' ')
disp('The Taylor series approximation of the nonlinear force vector')
fprintf('contains the following order(s): %s\n',num2str(nonlin_order));
disp(' ')
fprintf('f = \n');
disp(f); 
disp(' ')
disp('---------   Building up nonlinear coefficient matrices   --------- ')
disp(' ')

p_r = cell(1,max(nonlin_order));
for i=nonlin_order
    p_r{i} = combinator(numel(spv),i,'p','r');
end

vec = cell(numel(spv),1);
for i=1:numel(spv)
    vec{i} = repmat(i,t_num(i),1);
end

pos_f = cat(1,vec{:});
loc = zeros(numel(monomial_list),4);

for i=1:numel(monomial_list)
    id_dum = id_mon{i};
    loc(i,:) = [id_dum(1),pos_f(i),find(any(all(bsxfun(@eq,id_dum(2:end),p_r{id_dum(1)}),2),3)),coeff_list(i)];
end

for i=nonlin_order
    k = k+1;
    hwait=updatewaitbar(hwait,k,order);
    co = find(loc(:,1)==i);
    ii = loc(co,2);
    jj = loc(co,3);
    v = loc(co,4);
    Gcoef = sparse(ii,jj,v,numel(spv),numel(spv)^i);
    G{i} = T\Gcoef*nkron(i,T);
end

for i = iterate_remaining
    k = k+1;
    hwait=updatewaitbar(hwait,k,order);
    G{i} = 0;
end
    
close(hwait)
disp('Done')
varout = G;

end