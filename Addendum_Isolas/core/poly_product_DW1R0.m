function [H0,H0ij] = poly_product_DW1R0(W0,W0ij,R0_cell,R0ij_cell,numvar,corder,max_order,symbolic)  
% Comment: the 0th order can be computed with this function 
nvar = numvar;

if symbolic == 1 
    H_sub = sym(zeros([repmat(max_order+1,1,nvar),2])); 
else
    H_sub = zeros([repmat(max_order+1,1,nvar),2]); 
end


if corder == 0     
    for m = 1:nvar
        for j=1:2
            % Create unit index vector for each direction coordinate i 
            ei = ones(1,nvar);  
            ei(m) = ei(m)+1;
            ei_cell = sprintf('%i,',ei);
            ei_cell = ei_cell(1:end-1);

            ui = ones(1,nvar);  
            ui_cell = sprintf('%i,',ui);
            ui_cell = ui_cell(1:end-1);

            R0 = R0_cell{m};
            eval(strcat('H_sub(',ui_cell,',',num2str(j),') = H_sub(',ui_cell,',',num2str(j),')+W0(',ei_cell,',',num2str(j),')*R0(',ui_cell,');'))
        end
    end
else
    order = corder;  
    
    % Compute all coefficients for the current order
    ncoef = nch(order,nvar);
    i_dummy = [order,zeros(1,nvar-1)];
    index_corder = zeros(ncoef,nvar);
    index_corder(1,:) = i_dummy;

    for i=2:ncoef
        i_dummy = genlexd(i_dummy,nvar);
        index_corder(i,:) = i_dummy;
    end
    
    % Identify relevant terms in W from Wij for the current order
    W0ij_aorder = sum(W0ij,2);
    W0ij_corder = W0ij(W0ij_aorder<=order+1,:);

    for j=1:nvar

        R0 = R0_cell{j};
        R0ij = R0ij_cell{j};

        nonz_W0coef = [];
        nonz_R0coef = [];
        nonz_kcoef = [];

        ej = zeros(1,nvar);
        ej(j) = 1;
        index_corder_tilde = index_corder+repmat(ej,ncoef,1);

        for l=1:size(W0ij_corder,1) 
            
            % Compute k_tilde-m   (k_tilde = current order + ej, m = looping index)
            diff_coef = index_corder_tilde - repmat(W0ij_corder(l,:),ncoef,1);
            
            % Identify all terms R0(k_tilde-m) <> 0 
            [~,index_k,index_h] = intersect(diff_coef,R0ij,'rows','stable');

            % Build list of all nonzero terms for W0, R0 and the current order k
             if ~isempty(index_k)
                nonz_W0coef = [nonz_W0coef;repmat(W0ij_corder(l,:),numel(index_h),1)];
                nonz_R0coef = [nonz_R0coef;R0ij(index_h,:)];
                nonz_kcoef =  [nonz_kcoef; index_corder(index_k,:)];
             end
        end        
        
        for i=1:size(nonz_kcoef,1)

            i_min = j;
            k = nonz_kcoef(i,:);
            m = nonz_W0coef(i,:);
            k_m = nonz_R0coef(i,:);

            if  m(i_min)==0 

            else
                if sum(m<=k+ej)==nvar
                    
                    m_cell = sprintf('%i,',m+ones(size(m)));
                    m_cell = m_cell(1:end-1);
                    
                    k_m_cell = sprintf('%i,',k_m+ones(size(k_m)));
                    k_m_cell = k_m_cell(1:end-1);
                    
                    k_cell = sprintf('%i,',k+ones(size(k)));
                    k_cell = k_cell(1:end-1);
                    
                    for l=1:2
                        eval(strcat('Hadd = m(i_min)*W0(',m_cell,',',num2str(l),')*R0(',k_m_cell,');'));
                        eval(strcat('H_sub(',k_cell,',',num2str(l),') = H_sub(',k_cell,',',num2str(l),')+Hadd;'));
                    end                    
                end
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
    