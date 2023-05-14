function X = randisingMetropolis1(m,n,tau,gamma,maxiter);

% function X = randisingMetropolis1(m,n,tau,gamma,maxiter);
% or X = randisingMetropolis1(X0,tau,gamma,maxiter);
% m,n is the dimension of the Ising model
% X0 is an initial state, the outcome will have the same dimension as X0
% tau and gamma are the parameters of the model
% maxiter is the number of Gibbs steps.  (default is 100)
% Metropolis sampler

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

Xstep = X0;

maxii = 10;
maxjj = maxii;
changes = 0;
% proposal distribution generates new states on a subset of the random field.
% The subset has sides no larger than maxii x maxjj and is defined by its
% corners (i1,j1) and (i2,j2). These boundaries are random, uniformly
% distributed on the full random field
for step = 1:maxiter
   boundaries = ceil(rand(1,4).*[m m n n]);
   i1 = min(boundaries(1:2)); i2 = max(boundaries(1:2));
   i2 = min(i2,i1+maxii);
   j1 = min(boundaries(3:4)); j2 = max(boundaries(3:4));
   j2 = min(j2,j1+maxjj);
   Xprop = Xstep;
   Xprop(i1:i2,j1:j2) = (rand>0.5)*2-1;
   i1 = max(1,i1-1);
   i2 = min(m,i2+1);
   j1 = max(1,j1-1);
   j2 = min(n,j2+1);
   H2stepH = 0;
   if i1<i2,
      H2stepH = tau*(Xstep(i1:i2-1,j1:j2).*Xstep(i1+1:i2,j1:j2));
      H2stepH = sum(sum(H2stepH));
   end
   H2stepV = 0;
   if j1<j2,
      H2stepV = tau*(Xstep(i1:i2,j1:j2-1).*Xstep(i1:i2,j1+1:j2));
      H2stepV = sum(sum(H2stepV));
   end
   H2propH = 0;
   if i1<i2,
      H2propH = tau*(Xprop(i1:i2-1,j1:j2).*Xprop(i1+1:i2,j1:j2));
      H2propH = sum(sum(H2propH));
   end
   H2propV = 0;
   if j1<j2,
      H2propV = tau*(Xprop(i1:i2,j1:j2-1).*Xprop(i1:i2,j1+1:j2));
      H2propV = sum(sum(H2propV));
   end
   H2step = H2stepH + H2stepV;
   H2prop = H2propH + H2propV;

   H1step = gamma*sum(sum(Xstep(i1:i2,j1:j2)));
   H1prop = gamma*sum(sum(Xprop(i1:i2,j1:j2)));
   a = exp(-H1prop-H2prop)/exp(-H1step-H2step);
   A = rand;
   if A<a, Xstep = Xprop; changes = changes+1; end
end

changes
X = Xstep;
