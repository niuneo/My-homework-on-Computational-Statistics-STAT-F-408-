
function ans = randlaplace(m,n,l)

% randlaplace -- generates (pseudo) random numbers according to the Laplace or
%                double exponential distribution
%  Usage
%    ans = randlaplace(dim,l)
%    ans = randlaplace(m,n,l)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    l      intensity of Laplace distribution (parameter), see below
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    generates real (pos. or neg.) random numbers according to the distribution
%    pdf: f(t) = (l/2)*exp(-l*|t|)   (probability density function)
%  Note
%    randlaplace calls matlab's RAND function and therefore changes RAND's 
%    state.
%  Examples
%    b = randlaplace(2,4,3)
%    b = randlaplace(size(a),3)
%  See also
%    help randexp


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
ans = randexp(dim,l).*sign(rand(dim)-0.5);


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
