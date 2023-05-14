
function ans = randcauchy(m,n,a,b)

% randcauchy -- generates (pseudo) random numbers according to the Cauchy
%            distribution
%  Usage
%    ans = randcauchy(dim,a,b)
%    ans = randcauchy(m,n,a,b)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    a,b    mode and scale parameters in density/distribution (see below)
%           default values for a and b are 0 and 1 respectively.
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    generates real random numbers according to the distribution
%    pdf: f(x) = 1/(pi*b) * 1/(1+(x-a)^2/b^2)
%                (probability density function)
%    CDF: F(x) = 1/pi * atan((x-a)/b) + 1/2 (Cumulative Distribution Function)
%  Note
%    randcauchy calls matlab's RAND function and therefore changes RAND's
%    state.
%  Examples
%    X = randcauchy(2,4,3,2)
%    X = randcauchy(size(a),3,2)
%  See also

if (nargin == 0)
  dim = 1; a = 0; b = 1;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m; b = 1; a = 0;
  else
     dim = 1; a = m; b = 1;
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
x = b*tan(pi*(r-0.5))+a;
ans = x ;


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
