function varout = kronGK(n,m,K,G)

G_check = G{m};
G_check_bool = numel(G_check(:))>1;

if G_check_bool
    
    [~,x] = nsumk(m,n); 
    comb = x(all(x,2),:);

    S = zeros(size(K{1},1)^m,2^n); 
    
    for i = 1:size(comb,1)
       S = S + kronproduct(K{comb(i,:)});
    end

    varout = G{m}*S;
else
    varout = zeros(size(K{1},1),2^n);

end

end

