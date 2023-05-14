
function ans = randlogist(m,n,a,b)

% randlogist -- generates (pseudo) random numbers according to the logistic
%            distribution
%  Usage
%    ans = randlogist(dim,a,b)
%    ans = randlogist(m,n,a,b)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    a,b    parameters in density/distribution (see below)
%           default values for a and b are 0 and 1 respectively.
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    generates real random numbers according to the distribution
%    pdf: f(x) = (1/b) * exp(-(x-a)/b) / (1 + exp(-(x-a)/b)) 
%                (probability density function)
%    CDF: F(x) = 1 / (1 + exp(-(x-a)/b)) (Cumulative Distribution Function)
%    mean: E(X) = a
%    variance: var(X) = b^2*pi^2/3
%  Note
%    randlogist calls matlab's RAND function and therefore changes RAND's state.
%  Examples
%    X = randlogist(2,4,0,1)
%    X = randlogist(size(a),0,1)
%  See also

if (nargin == 0)
  dim = 1; a = 0; b = 1;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m; b = 1; a = 0;
  else
     dim = 1; a = 0; b = m;
  end
elseif (nargin == 2)
  if (length(m) > 1)
     dim = m; b = n; a = 0;
  else
     dim = [m n]; b = 1; a = 0;
  end
elseif (nargin == 3)
  if (length(m) > 1)
     dim = m; b = a; a = n;
  else
     dim = [m n]; b = a; a = 0;
  end
else
  dim = [m n];
end
r = rand(dim);
x = a + b * log(r./(1-r));
ans = x ;


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
