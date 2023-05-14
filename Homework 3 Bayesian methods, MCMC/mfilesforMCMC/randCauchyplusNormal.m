
function X = randCauchyplusNormal(m,n,mu,sigma,a,b)

% generates X ~ K*exp(-(x-mu)^2/(2*sigma^2))/(1+((x-a)/b)^2)


R = m*n;
X = zeros(1,R);
reject = (1:R);
while R > 0,
   Xr = randn(1,R)*sigma+mu;
   X(reject) = Xr;
   U = rand(size(Xr));
   fXoverMgX = 1./(1+(Xr-a).^2/b^2);
   reject = reject(U > fXoverMgX);
   R = length(reject);
end
X = reshape(X,m,n);
