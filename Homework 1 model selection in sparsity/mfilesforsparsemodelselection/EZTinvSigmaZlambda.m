
function s = EZTinvSigmaZlambda(Sigma,lambda,label,nsimul,nMCMCit);

% EZTinvSigmaZlambda -- MCMC sampling for 
%                       E(Z_I^T*Sigma^(-1)Z_I|A)
%                       with A the event abs(Z(i))>lambda)  <=> label(i)==1
%                            I = {i|label(i) = 1}
%                            Z ~ N(0,Sigma)
% Usage
%   s = EZTinvSigmaZlambda(Sigma,lambda,label,nsimul);
% Inputs
%   Sigma    square semi-positive definite matrix, correlation matrix
%   lambda   postive real number (scalar) (default 0)
%   label    binary vector  (default all 1's)
%   nsimul   number of simulations in MCMC (default 10000)
% Outputs
%   s        real number
% Description
%   This routines calls randmultivariatenormalabovethreshold for efficient
%   generation of normal-above-threshold random vectors
%   This routine is used in (among others) LARSell0
% See also
%   help randmultivariatenormalabovethreshold
%   help LARSell0

[n m] = size(Sigma);
if n~=m,
   error('Input covariance matrix must be square')
end
if max(max(abs(Sigma-Sigma')))>eps,
   error('Input covariance matrix must be symmetric')
end
[E Lambda] = eig(Sigma);
toll = 1.e-10;
if any(diag(Lambda)< -toll),
   error('Input covariance matrix must be semi-positive definite')
end

if nargin<2,
   lambda = 0;
end
if nargin<3,
   label = ones(n,1)
end
if nargin<4,
   nsimul = 10000;
end
if nargin<5,
   nMCMCit = 100;
end
if length(label) == 1,
   if label>1,
      nsimul = label;
      label = ones(n,1);
   end
end

I = find(abs(1-label)<toll);

if max(abs(1-label))<toll,
   V = randmultivariatenormalabovethreshold(nsimul,Sigma,lambda);
else
   V = randmultivariatenormalthreshold(nsimul,Sigma,lambda,label,nMCMCit);
   V = V(I,1:nsimul); Sigma = Sigma(I,I);
end
X = sum(V.*(Sigma\V),1);
Xbar = mean(X);
s = Xbar;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
