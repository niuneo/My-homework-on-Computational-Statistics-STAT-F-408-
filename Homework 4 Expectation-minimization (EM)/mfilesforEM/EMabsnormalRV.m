
function [muhat stdev2hat it] = EMabsnormalRV(X,varargin);

% EMabsnormalRV
%  Usage
%    [muhat stdev2hat it] = EMabsnormalRV(X);
%    [muhat stdev2hat it] = EMabsnormalRV(X,'maxit',maxit,'weights',w);
%  Inputs
%    X      (positive) observations; if X contains negative values, the
%                       absolute value is taken
%    'maxit',maxit  maximum number of iteration steps
%    'weights',w    weights of observations
%  Outputs
%    muhat
%    stdev2hat
%    it         number of iterations
%  This version uses expectation-maximization to find the expected values of
%  the signed observations
%  This version tries to reconstruct the original signed observations from the
%  absolute values

maxit = 100;
w = NaN;
nvarargin = length(varargin);
varargin1 = {varargin{:},NaN};
k = 1;
while k <= nvarargin,
   vark = varargin1{k};
   if ischar(vark), switch(vark)
   case {'maxit'}
      maxit = varargin{k+1};
      k = k+2;
   case {'weight','weights'}
      w = varargin{k+1};
      k = k+2;
   end
   elseif length(vark) == 1,
      maxit = vark;
      k = k+1;
   elseif length(vark) == length(X),
      w = vark;
      k = k+1;
   else
      k = k+1;
   end
end

X = column(abs(X));
if isnan(w), w = ones(size(X)); end

muhat = sum(w.*X)/sum(w);
% starting off with muhat = 0 does not work: for symmetry reasons, this is a
% local optimum of the log-likelihood curve
X2bar = sum(X.^2.*w)/sum(w);
stdev2hat = X2bar - muhat^2;

% expectation-maximization for variance estimator
it = 0;
muhat0 = Inf;
stdev2hat0 = 0;
while max(abs(stdev2hat-stdev2hat0)/stdev2hat,...
          abs(muhat-muhat0)/muhat)>1.e-6,
   it = it+1;
   % expectation
   stdevhat = sqrt(stdev2hat);
   % priorPneg = cumgauss(-muhat,stdevhat); priorPpos = 1-priorPneg;
   % postPneg = gauss(X+muhat,stdevhat)./...
   %           (gauss(X+muhat,stdevhat)+gauss(X-muhat,stdevhat));
   postPneg = exp(-2*abs(muhat)*X/stdev2hat)./...
             (exp(-2*abs(muhat)*X/stdev2hat)+1);
   if muhat < 0,
      postPpos = postPneg; postPneg = 1-postPpos;
   else
      postPpos = 1-postPneg;
   end
   % maximization (following line maximizes expected log-likelihood)
   muhat0 = muhat;
   stdev2hat0 = stdev2hat;
   muhat = sum(w.*postPpos.*X-w.*postPneg.*X)/sum(w);
   stdev2hat = sum(w.*postPneg.*(X+muhat).^2)/sum(w) + ...
               sum(w.*postPpos.*(X-muhat).^2)/sum(w);
   if it>maxit,
      muhat0 = muhat; stdev2hat0 = stdev2hat;
   end
end
