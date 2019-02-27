function varout = kronKR(n,m,K,R)

Rm = {eye(2),R{n+1-m}};
comb = eye(m)+ones(m);
S = zeros(2^m,2^n); 

for i = 1:size(comb,1)
   S = S + kronproduct(Rm{comb(i,:)});
end

varout = K{m}*S;

end

