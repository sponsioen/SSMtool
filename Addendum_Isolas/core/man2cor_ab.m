function  [Q,Qij] = man2cor_ab(W,T,cor_i,n_spv,n_var,max_order)    
    
    Q = zeros([repmat(max_order+1,1,n_var),2]);
    for i=1:n_spv
        Q = Q + T(cor_i,i)*W{i};
    end
    
    Qijdum =  cell(1,2);
    for i=1:2
        Qx1dum = Q(:,:,i);
        xi = find(Qx1dum);
        Xi_cell = cell(1,n_var);
        [Xi_cell{:}] = ind2sub(size(Qx1dum),xi);
        Xi = [Xi_cell{:}];
        Qijdum{i} =  (Xi-ones(size(Xi))).';
    end
    
    Qij = unique([Qijdum{:}].','rows');
    
end