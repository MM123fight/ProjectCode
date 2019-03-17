function SVMtoLP(experiment,k,lambda)
    [y,X] = libsvmread(char(experiment));
    [N,d] = size(X);
    m = (k-1)*N;
    n = 2*k*d+N;
    c = ones(n,1);
    c(1:2*k*d) = lambda*ones(2*k*d,1);
    bi = zeros(m,1);
    value = [];
    rowindex = [];
    colindex = [];
    for i = 1:N
        vi = X(y(i),1:d)';
        index_i = find(vi~=0);
        for j = 1:k
            if(j~=y(i))
                vj = X(j,1:d)';
                index_j = find(vj~=0);
                if(j < y(i))
                    value = [value;vj;-vj;-vi;vi;1 ];
                elseif(j>y(i))
                    value = [value;-vi;vi;vj;-vj;1 ];
                end
            end
        end
    end
    
end