function [H_CO,H_COij] = poly_power(U,Uij,Hp,Hpij,numvar,corder,power)
% Comment: the 0th order is incorporated into this function
nvar = numvar;
max_power = power;

% General recurrence formula: H_p(k) = (p/k_i)Sum_(m<=k, m_i>0)(m_i*U(m)H_(p-1)(k-m))
for p = 2:max_power
    
    % Create H_CO_sub 
    H_CO_sub = Hp{p};

    Hij = Hpij{p-1};
    H = Hp{p-1};

    % Compute 0th order term 
    ei = ones(1,nvar);  
    ei_cell = sprintf('%i,',ei);
    ei_cell = ei_cell(1:end-1);
    eval(strcat('H_CO_sub(',ei_cell,') = U(',ei_cell,')^p;'))

    for j = corder
        order = j; 
        
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
        Uij_corder = Uij(Uij_aorder<=j,:);
        
        nonz_Ucoef = [];
        nonz_Hcoef = [];
        nonz_kcoef = [];
        
        for l=1:size(Uij_corder,1) 
            
            % Compute k-m   (k = current order, m = looping index)
            diff_coef = index_corder - repmat(Uij_corder(l,:),ncoef,1);
            
            % Identify all terms H(k-m) <> 0 
            [~,index_k,index_h] = intersect(diff_coef,Hij,'rows');

            % Build list of all nonzero terms for U, V and the current order k
             if ~isempty(index_k)
                nonz_Ucoef = [nonz_Ucoef;repmat(Uij_corder(l,:),numel(index_h),1)];
                nonz_Hcoef = [nonz_Hcoef;Hij(index_h,:)];
                nonz_kcoef = [nonz_kcoef; index_corder(index_k,:)];
             end
        end

        % Loop over all found k's and add the result to H_sub(k)
        for i=1:size(nonz_kcoef,1)

            k = nonz_kcoef(i,:);
            k_i = min(k(k>0));
            
            i_min_dum = find(k==k_i);
            i_min = i_min_dum(1);

            m = nonz_Ucoef(i,:);
            k_m = nonz_Hcoef(i,:);

            if  m(i_min)==0

            else
                if sum(m<=k)==nvar
                    m_cell = sprintf('%i,',m+ones(size(m)));
                    m_cell = m_cell(1:end-1);
                    
                    k_m_cell = sprintf('%i,',k_m+ones(size(k_m)));
                    k_m_cell = k_m_cell(1:end-1);
                    
                    k_cell = sprintf('%i,',k+ones(size(k)));
                    k_cell = k_cell(1:end-1);
                    
                    eval(strcat('Hadd = (p/k_i)*m(i_min)*U(',m_cell,')*H(',k_m_cell,');'));
                    eval(strcat('H_CO_sub(',k_cell,') = H_CO_sub(',k_cell,')+Hadd;'));
                end
            end

        end
    end
    % Find nonzero elements in H_CO_sub
    xi = find(H_CO_sub);
    Xi_cell = cell(1,nvar);
    [Xi_cell{:}] = ind2sub(size(H_CO_sub),xi);
    Xi = [Xi_cell{:}];
    Hpij{p} = Xi-ones(size(Xi));
    Hp{p} = H_CO_sub;
end

% Create H_CO and H_COij updated with the current order 
H_COij = Hpij;
H_CO = Hp;
end




