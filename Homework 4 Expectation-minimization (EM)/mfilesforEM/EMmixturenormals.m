
function [p mu stdev pSpost] = EMmixturenormals(X,k,maxit,p,mu,stdev);

% EMmixturenormals
%  Usage
%    [p lambda mu stdev pSpost] = EMmixturenormals(X,k,maxit,p,mu,stdev);
%  Inputs
%    X      observations
%    k      number of hidden states
%    maxit  maximum number of iterations
%    p      vector of size k: initial guesses of P(S=s), s=1,...,k
%    mu     vector of size k: initial guesses of expected transformed RV
%    stdev  vector of size k: initial guesses of standard deviation of
%           transformed RV's
%  Outputs
%    p      vector of size k: estimate of P(S=s), s=1,...,k
%    mu     vector of size k: estimate of expected transformed RV
%    stdev  vector of size k: estimate of standard deviation of transformed RV
%    pSpost matrix of size [n k]: likelihood of hidden states for every 
%           observation

displayintermediate = true;
stopcrit = 1.e-5;

if nargin < 2, k = 2; end
if nargin < 3, maxit = 100; end
if nargin < 4, p = ones(1,k)/k; end
if nargin < 5,
   mu = zeros(1,k);
   n = length(X);
   q = [0 round(cumsum(p)*n)];
   Xsort = sort(X);
   for s=1:k,
     Xs = Xsort(q(s)+1:q(s+1));
     mu(s) = mean(Xs);
   end
end
if nargin < 6, 
   stdev = zeros(1,k);
   n = length(X);
   q = [0 round(cumsum(p)*n)];
   Xsort = sort(X);
   for s=1:k,
     Xs = Xsort(q(s)+1:q(s+1));
     stdev(s) = sqrt(mean(Xs.*Xs)-mu(s)^2);
   end
end



stopiter = false;
it = 1;
n = length(X);
X = column(X);
p = row(p);

mupast = mu;
while stopiter == false,
   stopiter = true;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % expectation: given global p, compute P(X_i=s)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   pSpost = postproblatentstatesmixturenormals(X,p,mu,stdev);
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % maximization: find estimators for p,mu,stdev
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   phat = sum(pSpost)/n;
   p = phat;
   for s = 1:k,
      ws = pSpost(1:n,s);
      [muhats stdev2hats] = EMabsnormalRV(X-min(X),ws);
      mu(s) = muhats+min(X);
      % muhats = sum(X.*ws)/sum(ws);
      % stdev2hats = sum((X-muhats).^2.*ws)/sum(ws);
      stdevhats = sqrt(stdev2hats);
      stdev(s) = stdevhats;
      if all(abs(mupast(1:it,s)-muhats)/stdevhats > stopcrit),
         stopiter = false;
      end
   end % for s = 1:k,
   mupast = [mupast; mu];

   if displayintermediate == true,
mupast
      p
      mu
      stdev
      it   
      disp('********************************************')
   end
   it = it+1;
   if it>maxit, stopiter = true; end
end 

% final estimation of posterior state probabilities
pSpost = postproblatentstatesmixturenormals(X,p,mu,stdev);

p = sum(pSpost)/n;
