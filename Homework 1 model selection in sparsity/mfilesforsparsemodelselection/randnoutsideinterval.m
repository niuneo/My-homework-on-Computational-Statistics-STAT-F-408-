
function ans = randnoutsideinterval(m,n,a,b);

% randnoutsideinterval -- generates (pseudo) random numbers 
%                         X = Z|Z outside [a,b] with Z ~ N(0,1)
%  Usage
%    ans = randnoutsideinterval(m,n,a,b)
%    ans = randnoutsideinterval(dim,a,b)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    a      left boundary of interval, default value is 0.
%    b      right boundary of interval, default value is -a.
%           (a and b are switched if a>b)
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    a and b are either scalars, or matrices of size equal to dim = [m n]
%  Note
%    randnoutsideinterval calls matlab's RANDN function and therefore changes 
%    RANDN's state.
%  Examples
%    X = randnoutsideinterval(2,4,1.4,5.3)
%    X = randnoutsideinterval(size(a),1.4,5.3)
%  See also
%    help invPhiUPhi
%    help randninsideinterval

if (nargin == 0)
  dim = 1;
  a = 0;
  b = 0;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m;
     a = 0;
     b = -a;
  else
     dim = [m 1];
     a = 0;
     b = -a;
  end
elseif (nargin == 2)
  if (length(m) > 1)
     dim = m;
     if length(n) == 2,
        a = n(1); b = n(2);
     elseif length(n) == 1,
        a = n;
        b = -a;
     else
        error('dimensions input incorrect')
     end
  else
     dim = [m n];
     a = 0;
     b = -a;
  end
elseif (nargin == 3)
  if (length(m) > 1)
     dim = m;
     b = a;
     a = n;
  else
     dim = [m n];
     if length(a) == 2,
        b = a(2); a = a(1);
     elseif length(a) == 1,
        b = -a;
     else
        error('dimensions input incorrect')
     end
  end
else
  dim = [m n];
end

if max(size(a)) > 1, if size(a) ~= dim, 
   error('size mismatch between a and dim = [m,n]')
end, end
if max(size(b)) > 1, if size(b) ~= dim,
   error('size mismatch between b and dim = [m,n]')
end, end

c = a;
s = (a>b);
a(s) = b(s);
b(s) = c(s);
% if a>b, ab = a; a = b; b = ab; end

% mn = dim(1)*dim(2);
if length(a) == 1, a = a*ones(dim); end
if length(b) == 1, b = b*ones(dim); end
poutside = cumgauss(a)+Rcumgauss(b);
o = find(poutside<eps); i = find(poutside>eps);
pleft = ones(size(a))/2;
pleft(i) = cumgauss(a(i))./poutside(i);
if ~isempty(o),
   aa = -a(o); bb = b(o);
   pleft(o) = exp(-aa.^2/2)./(exp(-aa.^2/2)+exp(-bb.^2/2)).*abs(bb./aa);
end
U = rand(dim);
left = find(U<=pleft);
right = find(U>pleft);
X = zeros(dim);
o = find(min(U,1-U).*poutside<eps);
X(left) = invcumgauss(U(left).*poutside(left));
X(right) = invcumgauss(1-(1-U(right)).*poutside(right));
if ~isempty(o),
   X(o) = invPhiUPhi(U(o),abs(a(o)),abs(b(o)));
   oleft = intersect(o,left);
   X(oleft) = - X(oleft);
end

ans = reshape(X,dim);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
