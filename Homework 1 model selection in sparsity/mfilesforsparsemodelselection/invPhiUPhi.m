
function x = invPhiUPhi(U,a,b,stdev);

% invPhiUPhi -- computes x = -invPhi(U*(Phi(-a)+Phi(-b))
%                        x = -invcumgauss(U*(cumgauss(-a)+cumgauss(-b))
%                        for large a and b (where cumgauss(-a) would give 0)
%                     or x = -invPhi(U*Phi(-a)+V*Phi(-b))
%  Usage
%    x = intPhiUPhi(U,a,b,stdev);
%    x = intPhiUPhi([U V],a,b,stdev);
%  Inputs
%    U      scalar, vector or matrix; values between 0 and 1
%    [U V]  scalars, vectors or matrices; values between 0 and 1
%    a      scalar or matrix of size(U)
%    b      scalar or matrix of size(U); default is Inf
%    stdev  default = 1
%  Outputs
%    x
%  Description
%  Note
%  Examples
%  See also
%    help cumgauss
%    help randnoutsideinterval 

if nargin < 4,
   stdev = 1;
end
if nargin < 3,
   b = Inf;
end
V = NaN;
if size(a,2) ~= size(U,2),
   ncolU = size(U,2); ncola = size(a,2);
   if ncolU == 2*ncola,
      nrowU = size(U,1);
      V = U(1:nrowU,ncola+1:ncolU);
      U = U(1:nrowU,1:ncola);
   end
end
if max(size(a))>1, if size(a) ~= size(U),
   error('size mismatch inputs U and a')
end, end
if max(size(b))>1, if size(b) ~= size(U),
   error('size mismatch inputs U and b')
end, end

if length(stdev)==1, stdev = stdev*ones(size(U)); end

if isnan(V),
   aa = min(a,b); bb = max(a,b); a = aa; b = bb;
   logc = log(1+(a./b).*exp((a.^2-b.^2)./(2*stdev.^2)));
else
   aa = min(a,b); bb = max(a,b); ii = find(aa~=a); a = aa; b = bb;
   UU = U; VV = V; VV(ii) = U(ii); UU(ii) = V(ii); U = UU; V = VV;
   logc = log(1+(V./U).*(a./b).*exp((a.^2-b.^2)./(2*stdev.^2)));
end
logU = log(U);
loga = log(a);
x1 = sqrt(a.^2+stdev.^2.*(loga - logU + logc));
x0 = zeros(size(x1));
i = find(abs(x0-x1)./abs(x1)>eps);
it = 0;
while length(i)>0,
   x0 = x1;
   x1 = sqrt(a.^2+stdev.^2.*(loga - log(x0) - logU + logc));
   it = it+1;
   i = find(abs(x0-x1)./abs(x1)>eps);
end
x = x1;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
