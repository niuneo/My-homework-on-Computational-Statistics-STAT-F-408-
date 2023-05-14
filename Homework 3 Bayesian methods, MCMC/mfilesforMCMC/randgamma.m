
function ans = randgamma(m,n,lambda,alfa)

% randgamma -- generates (pseudo) random numbers according to the gamma 
%               distribution
%  Usage
%    ans = randgamma(m,n,lambda,alfa)
%    ans = randgamma(dim,lambda,alfa)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    lambda parameter of Gamma distribution, default is 1. Lambda is the
%           intensity of the underlying Poisson process (i.e., the expected
%           number of arrivals per time unit)
%    alfa   parameter of Gamma distribution, default value is 1. 
%           If r=1, then the Gamma distribution coincides with the exponential
%           distribution. 
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    generates real random numbers according to the distribution
%    pdf: f(x) = lambda^alfa*exp(-lambda*x)*x^(alfa-1)/Gamma(alfa)
%    expected value: E(X) = alfa/lambda
%    variance: var(X) = alfa/lambda^2
%  Note
%    randgamma calls matlab's RAND function and therefore changes RAND's
%    state.
%  Examples
%    X = randgamma(2,4,1.4,5.3)
%    X = randgamma(size(a),1.4,5.3)
%  See also
%    help randexp
%    help randerlang

if (nargin == 0)
  dim = 1; lambda = 1; alfa = 1;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m; lambda = 1; alfa = 1;
  else
     dim = 1; alfa = 1; lambda = m;
  end
elseif (nargin == 2)
  if (length(m) > 1)
     dim = m; lambda = n; alfa = 1;
  else
     dim = [m n]; lambda = 1; alfa = 1;
  end
elseif (nargin == 3)
  if (length(m) > 1)
     dim = m; alfa = lambda; lambda = n;
  else
     dim = [m n]; alfa = 1;
  end
else
  dim = [m n];
end

m = dim(1); n = dim(2);

r = floor(alfa);
ralfa = alfa - r;
y = randexp([m n r],lambda);
x = sum(y,3);

if ralfa > eps,
   e = exp(1);
   p = e/(e+ralfa);

   mn = m*n;

   xi = zeros(1,mn);
   k = (1:mn);
   K = mn;

   while K > 0,
      % We generate X ~ gX(x)
      %   with gX(x) = p*ralfa*x^(ralfa-1) for 0<x<1
      %   and  gX(x) = (1-p)*exp(1-x) for x>1
      % X = (V<p) * U*(1/ralfa) +
      %     (V>p) * (1-log(U))
      % We accept X if fX(X)/[M*gX(X)] < R, with R ~ uniform[0,1]
      %
      V = rand(1,K);
      U = rand(1,K);
      R = rand(1,K);
      I = find(V<=p); II = find(V>p);
      X = zeros(1,K);
      X(I) = U(I).^(1/ralfa); X(II) = 1-log(U(II));
      fXoverMgX = zeros(1,K);
      fXoverMgX(I) = exp(-X(I));
      fXoverMgX(II) = X(II).^(ralfa-1);
      % gX = zeros(1,K);
      % gX(I) = p*ralfa*X(I).^(ralfa-1);
      % gX(II) = (1-p)*exp(1-X(II));
      reject = find(R>fXoverMgX);
      xi(k) = X;
      k = k(reject);
      K = length(k);
   end
   xi = reshape(xi,m,n)/lambda;
else
   xi = zeros(size(x));
end

ans = x+xi ;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
