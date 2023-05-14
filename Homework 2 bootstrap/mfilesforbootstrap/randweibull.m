
function ans = randweibull(m,n,l,r)

% randweibull -- generates (pseudo) random numbers according to a Weibull law
%  Usage
%    ans = randweibull(dim,l,r)
%    ans = randweibull(m,n,l,r)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    l,r    parameters of the Weibull distribution, see below
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    generates positive real random numbers according to the distribution
%    pdf: f(t) = r*l*(l*t)^(r-1)*exp((-l*t)^r) (probability density function)
%    CDF: F(t) = 1 - exp(-(l*t)^r) (Cumulative Distribution Function)
%  Note
%    randweibull calls matlab's RAND function and therefore changes RAND's state. 
%  Examples
%    b = randweibull(2,4,3,2)
%    b = randweibull(size(a),3,2)
%  See also
%    help randexp

if (nargin == 0)
  dim = 1; l = 1; r = 1;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m; l = 1;
  else
     dim = 1; l = m;
  end
  r = 1;
elseif (nargin == 2)
  if (length(m) > 1)
     dim = m; l = n; r = 1;
  else
     dim = 1; l = m; r = n;
  end
elseif (nargin == 3)
  if (length(m) > 1)
     dim = m; r = l; l = n;
  else
     dim = [m n]; r = 1; 
  end
else
  dim = [m n];
end
U = rand(dim);
x = abs(log(U)).^(1/r)/l;
ans = x ;


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
