
function ans = randt(m,n,dof)

% randt -- generates (pseudo) random numbers according to the Student's t
%          distribution
%  Usage
%    ans = randt(dim,dof)
%    ans = randt(m,n,dof)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    dof    number of degrees of freeddom, parameter of t distribution
%           default values is 1.
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    ans is generated as 
%    ans = randn(dim)./sqrt(randchi2(dim,dof)/dof);
%    where randn(dim) and randchi2(dim,dof) are independent random vectors of
%    equal size
%  Note
%    randt calls matlab's RANDN function and therefore changes RANDN's
%    state.
%  Examples
%    X = randt(2,4,3)
%    X = randt(size(a),3)
%  See also
%    help randchi2

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

ans = randn(dim)./sqrt(randchi2(dim,dof)/dof);


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
