
function [betahat, muhat, x, Cp, lambda, Aopt] = LARS(Y,K,stdev,n1)

% LARS - Least Angle Regression for LASSO (beta-version; subject to change)
%  Usage
%    [betahat muhat x Cp lambda Aopt] = LARS(Y,K,stdev,n1)
%  Inputs
%    Y      observations in sparsity model Y = K*beta+noise
%    K      model of covariates (rectangular matrix)
%    stdev  (estimated) standard deviation (default is 1; current version does
%           not estimate stdev)
%    n1     desired number of nonzeros (default is n=length(Y))
%  Outputs
%    betahat estimator beta
%    muhat   estimator mu = K*beta
%    x       ordered set of selected variables (indices in {1,..,length(beta)})
%    Cp      Mallows' Cp values of nested selections
%    lambda  values of penalty parameter corresponding to x and Cp
%    Aopt    set of selected variables leading to minimum Cp value
%  Description
%  Note
%     see p.74 in Elements of statistical learning, second edition;
%     downloadable from http://www-stat.stanford.edu/~tibs/ElemStatLearn/
%     see also paper LARS, Ann.Stat.2004, pages 407-499
%  See also:
%    help iterativeST

if nargin < 4, n1 = n; end
if nargin < 3, stdev = 1; end

Y = column(Y);

ybar = mean(Y);
Y = Y-ybar;
muhat = zeros(size(Y));
n = length(Y);
m = size(K,2);
colnormK = sqrt(sum(K.^2,1));
X = K*diag(1./colnormK);

k = 0;
dof = k;
% Cp = [Cp(-1) Cp(0) Cp(1)]
Cp = [Inf, norm(Y-muhat)^2/n+2*dof/n*stdev^2-stdev^2];
Cpmin = Inf;
muhatopt = zeros(size(Y)); Aopt = [];
lambda = [Inf, Inf];
A = []; x = A;
Aprime = (1:m);
% while Cp(k+1) > Cp(k+2),
while k < n1,
   muhat0 = muhat;
   chat = X'*(Y-muhat0);
   [Chat, j] = max(abs(chat));
   if k>0,
      sA = sign(chat(A));
      XA = X(1:n,A)*diag(sA);
      GA = XA'*XA;
      AA = 1/sqrt(ones(1,k)*(GA\ones(k,1)));
      wA = AA*(GA\ones(k,1));
      uA = XA*wA;
   else
      uA = X(1:n,j); uA = uA/norm(uA);
      AA = 1;
   end
   a = X'*uA;
   gamma = [(Chat-chat(Aprime))./(AA-a(Aprime)), ...
            (Chat+chat(Aprime))./(AA+a(Aprime))];
   gamma(gamma<0) = Inf;
   gamma = min(gamma,[],2);
   [gammahat, j] = min(gamma);
   muhat = muhat0 + gammahat*uA;
   newAprime = [Aprime(1:j-1), Aprime(j+1:m-k)]; 
   j = Aprime(j); A = [A, j]; Aprime = newAprime;
   k = k+1;
   dof = k;
   lambda = [lambda, Chat];
   Chat = Chat -gammahat*AA;
   Cpk = norm(Y-muhat)^2/n+2*dof/n*stdev^2-stdev^2;
   if Cpk < Cpmin,
      Cpmin = Cpk;
      Aopt = A;
      muhatopt = muhat;
   end
   Cp = [Cp, Cpk];
end
muhat = muhatopt+ybar;
A = Aopt;
XA = X(1:n,A);

% HT
% betahatA = (XA'*XA)\(XA'*Y);
% ST
betahatA = (XA'*XA)\(XA'*muhat);
betahat = zeros(1,m);
betahat(A) = betahatA;

x = A;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
