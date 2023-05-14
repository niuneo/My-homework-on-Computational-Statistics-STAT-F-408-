
function pSpost = postproblatentstatesmixturenormals(x,pSprior,mu,stdev);

% postproblatentstatesmixturenormals.m
%  finds P(S=s|X=x) in a model where S ~ pSprior and
%        X|S=s ~ N(mu(s),stdev(s))
%  Usage
%    pSpost = postproblatentstatesmixturenormals(x,pSprior,mu,stdev);
%  Inputs
%    X       vector of observations
%    pSprior vector of size k with sum(pSprior)=1 and positive entries; 
%            prior model for S
%    mu      vector of size k; parameters for conditional models X|S=s
%    stdev   vector of size k; parameters for conditional models X|S=s
%  Outputs
%    pSpost  posterior probabilities P(S=s|X=x) (vector of size(X))

k = length(pSprior);
n = length(x);
pSprior = row(pSprior); x = column(x);

pSpost = ones(n,k)/k;

if k>1,
   pSpost = zeros(n,k);
   for s = 1:k,
      fXgivenS = gauss(x-mu(s),stdev(s));
      pSpost(1:n,s) = pSprior(s)*fXgivenS;
   end
   fX = sum(pSpost')';
   pSpost = pSpost./repmat(fX,1,k);
end
