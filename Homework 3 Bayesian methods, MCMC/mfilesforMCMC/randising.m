
function X = randising(m,n,tau,gamma,maxiter);

% function X = randising(m,n,tau,gamma,maxiter);
% or X = randising(X0,tau,gamma,maxiter);
%
% GIBBS sampler for 2-dim. Ising model
%
% m,n is the dimension of the Ising model
% X0 is an initial state, the outcome will have the same dimension as X0
% tau and gamma are the parameters of the model
% maxiter is the number of Gibbs steps.  (default is 100)
% MCMC with quincunx scheme

if nargin < 5, maxiter = 100; end
if nargin < 4, gamma = 0; end
if nargin < 3, tau = 0; end

if any(size(m) ~= [1 1]),
   if nargin > 3, maxiter = gamma; end
   if nargin > 2, gamma = tau; end
   if nargin > 1, tau = n; end
   if all(size(m) == [1 2]) | all(size(m) == [2 1]),
      n = m(2);
      m = m(1);
      U = rand(m,n);
      X0 = (U<1/2)*2-1;
   else
      X0 = m;
      m = size(X0,1); n = size(X0,2);
   end
else
   U = rand(m,n);
   X0 = (U<1/2)*2-1;
end
drift1 = exp(-gamma)/(exp(-gamma)+exp(gamma));
drift0 = 1-drift1;

X = X0;
X0 = [zeros(1,n+2); zeros(m,1) X zeros(m,1); zeros(1,n+2);];

j = ones(m,1)*(1:n);
i = (1:m)'*ones(1,n);
quincunx1 = i+j - floor((i+j)/2)*2;
quincunx2 = ones(size(quincunx1))-quincunx1;

for iter = 1:maxiter
   Uiter = rand(m,n);
   H1 = tau*(X0(1:m,2:n+1)+X0(3:m+2,2:n+1)+X0(2:m+1,1:n)+X0(2:m+1,3:n+2));
   p1iter = drift1*exp(-H1)./(drift0*exp(H1)+drift1*exp(-H1));
   Xiter1 = (Uiter < p1iter)*2-1;
   Xiter = X0(2:m+1,2:n+1);
   Xiter = quincunx1.*Xiter1+quincunx2.*Xiter;
   X0(2:m+1,2:n+1) = Xiter;
   H1 = tau*(X0(1:m,2:n+1)+X0(3:m+2,2:n+1)+X0(2:m+1,1:n)+X0(2:m+1,3:n+2));
   p1iter = drift1*exp(-H1)./(drift0*exp(H1)+drift1*exp(-H1));
   Xiter2 = (Uiter < p1iter)*2-1;
   Xiter = quincunx2.*Xiter2+quincunx1.*Xiter;
   X0(2:m+1,2:n+1) = Xiter;
end
   
X = Xiter;
