
function ans = randchi2(m,n,dof)

% randchi2 -- generates (pseudo) random numbers according to the chi-square
%             distribution
%  Usage
%    ans = randchi2(dim,dof)
%    ans = randchi2(m,n,dof)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    dof    number of degrees of freeddom, parameter of chi-square distribution
%           default values is 1.
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%  Note
%    randchi2 calls matlab's RANDN function and therefore changes RANDN's
%    state.
%  Examples
%    X = randchi2(2,4,3)
%    X = randchi2(size(a),3)
%  See also
%    help randt

if (nargin == 0)
  dim = 1; dof = 1;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m; dof = 1;
  else
     dim = 1; dof = m;
  end
elseif (nargin == 2)
  if (length(m) > 1)
     dim = m; dof = n;
  else
     dim = [m n]; dof = 1;
  end
else
     dim = [m n];
end
x = zeros(dim);
for i = 1:dof,
   x = x + randn(dim).^2;
end
ans = x ;


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
