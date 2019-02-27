function varout = kronGK1n(n,K,G)
    varout = G{n}*nkron(n,K{1});
end

