
function betahat = iterativeST(Y,K,thr,maxit,betahat0);
% vhats = iterativeST(w,K,thr(t),600)
% iterativeST -- Iterative soft thresholding
%  Usage
%    betahat = iterativeST(Y,K,thr,maxit,betahat0);
%  Inputs
%    Y      observations in sparsity model Y = K*beta+noise 
%    K      model of covariates (rectangular matrix)
%    thr    fixed threshold used in iterations
%    maxit  maximum number of iterations (default is 1000)
%    betahat0  initial value for betahat (default is zero)
%  Outputs
%    betahat   estimate of sparse vector beta
%  Description
%    Iterative soft thresholding:
%    betahat_new = ST(betahat_old + K'*(Y-K*betahat_old),thr)
%  Note
%  See also:
%    help iterativeHT
%    help iterativeGCVST
%    help iterativeMSEST
%    help LARS

ysize = size(Y);
Y = column(Y);

if nargin < 2,
   K = eye(length(Y));
end
[n m] = size(K);
if nargin < 3,
   thr = sqrt(2*log(m));
end

if nargin < 4,
   maxit = 1000;
  
end

if nargin < 5,
   betahat0 = zeros(m,1);
end

beta0 = ones(m,1);
beta1 = betahat0;

if norm(beta1-beta0) < eps, % i.e., if betahat0 happens to be ones(m,1)
   beta0 = zeros(m,1);
end

lambda = thr;

i = 1;
while norm(beta1-beta0) > eps,
   beta0 = beta1;
   y1 = beta0 + K'*(Y-K*beta0);
   beta1 = ST(y1,lambda);
%    ST(x,t) = (abs(x)-t).*sign(x).*(abs(x)>t)
%            = x.*(abs(x)>t).*(1-t./abs(x));
   i = i+1; 
   if i>maxit, 
      beta0 = beta1; 
      warning('iterativeST: maximum number of iterations reached')
   end
end
betahat = beta1;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 

