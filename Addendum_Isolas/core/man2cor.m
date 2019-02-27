function  [P,Pij] = man2cor(W,T,cor_i,n_spv,n_var,max_order)

    P = zeros(repmat(max_order+1,1,n_var));
    for i=1:n_spv
        P = P + T(cor_i,i)*W{i};
    end
    
    xi = find(P);
    Xi_cell = cell(1,n_var);
    [Xi_cell{:}] = ind2sub(size(P),xi);
    Xi = [Xi_cell{:}];
    Pij =  Xi-ones(size(Xi));

end