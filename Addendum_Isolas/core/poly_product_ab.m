function [H0,H0ij] = poly_product_ab(U,Uij,V,Vij,H,numvar,corder)
% Comment: the 0th order can be computed with this function 
nvar = numvar;
H_sub = H;
   
for n = corder
    
    order = n;  
    
    % Compute all coefficients for the current order
    ncoef = nch(order,nvar);
    i_dummy = [order,zeros(1,nvar-1)];
    index_corder = zeros(ncoef,nvar);
    index_corder(1,:) = i_dummy;
    
    for i=2:ncoef
        i_dummy = genlexd(i_dummy,nvar);
        index_corder(i,:) = i_dummy;
    end
    
    % Identify relevant terms in U from Uij for the current order
    Uij_aorder = sum(Uij,2);
    Uij_corder = Uij(Uij_aorder<=order,:);

    nonz_Ucoef = [];
    nonz_Vcoef = [];
    nonz_kcoef = [];

    for l=1:size(Uij_corder,1) 
        
        % Compute k-m   (k = current order, m = looping index)
        diff_coef = index_corder - repmat(Uij_corder(l,:),ncoef,1);
        
        % Identify all terms V(k-m) <> 0 
        [~,index_k,index_h] = intersect(diff_coef,Vij,'rows','stable');

        % Build list of all nonzero terms for U, V and the current order k
        if ~isempty(index_k)
            nonz_Ucoef = [nonz_Ucoef;repmat(Uij_corder(l,:),numel(index_h),1)];
            nonz_Vcoef = [nonz_Vcoef;Vij(index_h,:)];
            nonz_kcoef = [nonz_kcoef; index_corder(index_k,:)];
        end
    end

    % Loop over all found k's and add the result to H_sub(k)
    for i=1:size(nonz_kcoef,1)

        k = nonz_kcoef(i,:);
        m = nonz_Ucoef(i,:);
        k_m = nonz_Vcoef(i,:);

        if sum(m<=k)==nvar
            
            m_cell = sprintf('%i,',m+ones(size(m)));
            m_cell = m_cell(1:end-1);
            
            k_m_cell = sprintf('%i,',k_m+ones(size(k_m)));
            k_m_cell = k_m_cell(1:end-1);
            
            k_cell = sprintf('%i,',k+ones(size(k)));
            k_cell = k_cell(1:end-1);
            
            for l=1:2
                eval(strcat('Hadd = U(',m_cell,')*V(',k_m_cell,',',num2str(l),');'));
                eval(strcat('H_sub(',k_cell,',',num2str(l),') = H_sub(',k_cell,',',num2str(l),')+Hadd;'));
            end
        end

    end

end


% Find nonzero elements in H_sub and list them in H0ij
H0ijdum = cell(1,2);

for i=1:2
    H_subdum = H_sub(:,:,i);
    xi = find(H_subdum);
    Xi_cell = cell(1,nvar);
    [Xi_cell{:}] = ind2sub(size(H_subdum),xi);
    Xi = [Xi_cell{:}];
    H0ijdum{i} =  (Xi-ones(size(Xi))).';
end

H0ij = unique([H0ijdum{:}].','rows');

% Create H0 updated with the current order 
H0 = H_sub;

end
    

