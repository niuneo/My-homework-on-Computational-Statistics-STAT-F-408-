n = 10000;
p = [0.2 0.7 0.1];
mu = [0.01 4 6];
stdev = [1 0.4 1];
k = length(p);

U = rand(n,1); pp = cumsum(p); S = ones(n,1);
for s = 2:k, S = S + (U>pp(s-1)); end

Z = randn(n,1);
W = zeros(n,1);
X = zeros(n,1);
for s = 1:k,
   ii = find(S==s);
   X(ii) = mu(s)+stdev(s)*Z(ii);
end

t = (0:k-1)/k;
rgb = [1 0 0; 0 1 0; 0 0 1;1 0 0];

maxit = 15;
[phat muhat stdevhat pSpost] = EMmixturenormals(X,3,maxit);


