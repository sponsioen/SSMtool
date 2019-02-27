function [vargout] = genlexd(c,k)
% Compute the next vector in descending lexicographical order. 
is = c(k);

for nv = k-1:-1:1
    is = is + c(nv);
    
    if c(nv)~=0
        c(nv)=c(nv)-1;
        c(nv+1)=is-c(nv);
        vargout = c;
        return;
    else
        c(nv+1)=0;
    end
end

c(1)=is;
vargout = c;


