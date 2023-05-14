
function [LL stdev2hat it] = ...
          profileloglikelihoodabsnormalfixedpoint(mu,X,maxit);

% profileloglikelihoodabsnormalfixedpoint
%  Usage
%    LL = profileloglikelihoodabsnormalfixedpoint(mu,X,maxit);
%    [LL stdev2hat it] = profileloglikelihoodabsnormalfixedpoint(mu,X,maxit);
%  Inputs
%    mu     (vector of) expected values of a normal random variable
%    X      observations; if X contains no negative values, they are considered
%                         as absolute values of normal random variables.
%                         Otherwise they are consider as proper normal rv's.
%    maxit  maximum number of iterations in fixed point optimization method
%  Outputs
%    LL         profile likelihood of mu
%    stdev2hat  variance estimator for given mu
%    it         number of iterations
% This version uses a fixed point method (i.e., numerical iteration) to
% optimize the loglikelihood (more precisely to find zeros of the partial
% derivative of the LL to stdev^2)
% In case X contains the observations with signs, we use them for initial
% estimation, but then run the iterations for absoute values

if nargin<3, maxit = 100; end
X = column(X);
mu = row(mu);
n = length(X);
m = length(mu);
X = repmat(X,1,m);
mu = repmat(mu,n,1);

if min(X(1:n,1))<-eps,
   stdev2hat = mean((X-mu).^2,1); stdev2hat = repmat(stdev2hat,n,1);
   X = abs(X);
else
   X2bar = mean(X(1:n,1).^2);
   stdev2hat = max(0,X2bar - mu.^2);
end
% Now follows the iteration: even if we know sign(X), we pretend not so, in
% order to check whether the optimum for sign(X) known is invariant under
% the iteration
% stdev2hat = mean((X-mu).^2,1); stdev2hat = repmat(stdev2hat,n,1);
it = 0;
stdev2hat1 = stdev2hat(1,1:m);
stdev2hat0 = Inf*ones(size(stdev2hat1));
while max(abs(stdev2hat1-stdev2hat0))/min(stdev2hat1)>1.e-6,
   it = it+1;
   stdevhat = sqrt(stdev2hat);
   D = gauss(X-mu,stdevhat)+gauss(X+mu,stdevhat);
   N = gauss(X-mu,stdevhat).*(X-mu).^2 + ...
       gauss(X+mu,stdevhat).*(X+mu).^2;
   stdev2hat = sum(N./D,1)/n;
   stdev2hat = repmat(stdev2hat,n,1);
   stdev2hat0 = stdev2hat1;
   stdev2hat1 = stdev2hat(1,1:m);
   if it>maxit,
      stdev2hat0 = stdev2hat1;
   end
end
stdevhat = sqrt(stdev2hat);

if min(X(1:n,1))>-eps,
   LL = log(gauss(X-mu,stdevhat)+gauss(X+mu,stdevhat));
else
   LL = log(gauss(X-mu,stdevhat));
end
LL = sum(LL,1);
stdev2hat = stdev2hat(1,1:m);
