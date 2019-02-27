function varout = kronproduct(varargin)

K = kron(varargin{1},varargin{2});


for i = 3:length(varargin)
    K = kron(K, varargin{i});
end

varout = K;

end




