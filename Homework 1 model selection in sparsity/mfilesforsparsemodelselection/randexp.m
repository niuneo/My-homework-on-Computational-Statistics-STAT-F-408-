
function ans = randexp(m,n,l)

% randexp -- generates (pseudo) random numbers according to exponential
%            distribution
%  Usage
%    ans = randexp(dim,l)
%    ans = randexp(m,n,l)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    l      intensity of exponential distribution (parameter), see below
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    generates positive real random numbers according to the distribution
%    pdf: f(t) = l*exp(-l*t)   (probability density function)
%    CDF: F(t) = 1 - exp(-l*t) (Cumulative Distribution Function)
%  Note
%    randexp calls matlab's RAND function and therefore changes RAND's state. 
%  Examples
%    b = randexp(2,4,3)
%    b = randexp(size(a),3)
%  See also
%    help exptimes for exponentially distributed data with varying intensities

if (nargin == 0)
  dim = 1; l = 1;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m; l = 1;
  else
     dim = 1; l = m;
  end
elseif (nargin == 2)
  if (length(m) > 1)
     dim = m; l = n;
  else
     dim = [m n]; l = 1;
  end
else
     dim = [m n];
end
U = rand(dim);
x = abs(log(U))/l;
ans = x ;


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
