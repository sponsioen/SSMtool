function varout = nkron(n,K)

Kr = kron(K,K);

for i = 3:n
    Kr = kron(Kr, K);
end

varout = Kr;

end




