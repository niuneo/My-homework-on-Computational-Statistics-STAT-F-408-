
function ans = randgumbel(m,n,a,b)

% randgumbel -- generates (pseudo) random numbers according to the Gumbel
%            distribution
%  Usage
%    ans = randgumbel(dim,a,b)
%    ans = randgumbel(m,n,a,b)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    a,b    parameters in density/distribution (see below)
%           default values for a and b are 1 and 0 respectively.
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    generates real random numbers according to the distribution
%    CDF: F(x) = exp(-exp(-a(x-b))) (Cumulative Distribution Function)
%  Note
%    randgumbel calls matlab's RAND function and therefore changes RAND's
%    state.
%  Examples
%    X = randgumbel(2,4,1,0)
%    X = randgumbel(size(a),1,0)
%  See also

if (nargin == 0)
  dim = 1; a = 1; b = 0;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m; b = 0; a = 1;
  else
     dim = 1; a = 1; b = m;
  end
elseif (nargin == 2)
  if (length(m) > 1)
     dim = m; b = n; a = 1;
  else
     dim = [m n]; b = 0; a = 1;
  end
elseif (nargin == 3)
  if (length(m) > 1)
     dim = m; b = a; a = n;
  else
     dim = [m n]; b = a; a = 1;
  end
else
  dim = [m n];
end
r = rand(dim);
x = b - log(abs(log(r)))/a;
ans = x ;


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
