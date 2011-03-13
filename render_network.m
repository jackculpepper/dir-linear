function array = render_network(A, m, clim)

[L M] = size(A);

n = M/m;

sz = sqrt(L);

buf = 1;

array = -ones(buf+m*(sz+buf), buf+n*(sz+buf));

k = 1;

for i = 1:m
    for j = 1:n
        array(buf+(i-1)*(sz+buf)+[1:sz], buf+(j-1)*(sz+buf)+[1:sz]) = ...
            reshape(A(:,k),sz,sz)/clim;
        k = k+1;
    end
end

