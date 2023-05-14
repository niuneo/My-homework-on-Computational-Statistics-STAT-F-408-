
function ans = randbeta(m,n,a,b);

% randbeta -- generates (pseudo) random numbers according to the beta 
%               distribution
%  Usage
%    ans = randbeta(m,n,a,b)
%    ans = randbeta(dim,a,b)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    a      parameter of Beta distribution, default value is 1. 
%    b      parameter of Beta distribution, default value is a.
%           distribution. 
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%  Note
%    randbeta calls matlab's RAND function and therefore changes RAND's
%    state.
%  Examples
%    X = randbeta(2,4,1.4,5.3)
%    X = randbeta(size(a),1.4,5.3)
%  See also
%    help randgamma

if (nargin == 0)
  dim = 1; a = 1; b = 1;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m; a = 1; b = 1;
  else
     dim = [m 1]; a = 1; b = 1;
  end
elseif (nargin == 2)
  if (length(m) > 1)
     dim = m; a = n; b = a;
  else
     dim = [m n]; a = 1; b = 1;
  end
elseif (nargin == 3)
  if (length(m) > 1)
     dim = m; b = a; a = n;
  else
     dim = [m n]; b = a;
  end
else
  dim = [m n];
end

X = randgamma(dim,1,a);
Y = randgamma(dim,1,b);
ans = X./(X+Y);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
