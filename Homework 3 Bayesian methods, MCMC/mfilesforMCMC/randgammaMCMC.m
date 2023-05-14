function ans = randgammaMCMC(m,n,lambda,alfa)

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
%    The values of X are generated as a Markov Chain, using a
%    Metropolis-Hastings sampler
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
mn = m*n;
X = zeros(mn,1);

r = floor(alfa);
ralfa = alfa - r;
Y = randexp(1,r,lambda);  % why take r samples?
X(1) = sum(Y);  

for k=2:mn,
   reject = true;
   while reject,
      Y = randexp(1,r,lambda); X(k) = sum(Y);  % explain?
      V = (X(k)/X(k-1))^ralfa;
      U = rand;
      if U<V, reject = false; end
   end
end
ans = reshape(X,m,n);
