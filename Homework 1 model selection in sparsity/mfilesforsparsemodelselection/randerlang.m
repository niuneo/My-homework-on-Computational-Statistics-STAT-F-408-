
function ans = randerlang(m,n,lambda,r)

% randerlang -- generates (pseudo) random numbers according to the gamma 
%               distribution with integer parameter (Erlang distribution)
%  Usage
%    ans = randerlang(m,n,lambda,r)
%    ans = randerlang(dim,lambda,r)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    lambda parameter of Erlang distribution, default is 1. Lambda is the
%           intensity of the underlying Poisson process (i.e., the expected
%           number of arrivals per time unit)
%    r      parameter of Erlang distribution, default value is 1. r is the
%           number of arrivals waited for. If r=1, then the Erlang distribution
%           coincides with the exponential distribution. r must be integer. Use
%           randgamma for non-integer r.
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    X ~ Erlang(lambda,r) <=> X equals (in distribution) the sum of r
%    independent, identically distributed, exponential random variables X_i, 
%    whose cumulative distribution functions are 1-exp(lambda*x)
%    E(X) = r/lambda  and  var(X) = r/lambda^2
%  Note
%    randerlang calls matlab's RAND function and therefore changes RAND's
%    state.
%  Examples
%    X = randerlang(2,4,1.4,5)
%    X = randerlang(size(a),1.4,5)
%  See also
%    help randexp
%    help randgamma

if (nargin == 0)
  dim = 1;
  lambda = 1;
  r = 1;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m;
     lambda = 1;
     r = 1;
  else
     dim = 1;
     r = 1;
     lambda = m;
  end
elseif (nargin == 2)
  if (length(m) > 1)
     dim = m;
     lambda = n;
     r = 1;
  else
     dim = [m n];
     lambda = 1;
     r = 1;
  end
elseif (nargin == 3)
  if (length(m) > 1)
     dim = m;
     r = lambda;
     lambda = n;
  else
     dim = [m n];
     r = 1;
  end
else
  dim = [m n];
end
if abs(r-round(r))>eps,
  error('Erlang not defined for non-integer parameters. Use randgamma instead')
end

m = dim(1); n = dim(2);

y = randexp([m n r],lambda);
x = sum(y,3);
ans = x ;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
