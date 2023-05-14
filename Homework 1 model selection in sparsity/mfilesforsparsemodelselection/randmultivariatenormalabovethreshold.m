
function ans = randmultivariatenormalabovethreshold(m,Sigma,thr,maxit)

% randmultivariatenormalabovethreshold -- generates (pseudo) random vectors 
%                         Y = X|all abs(X_i) > thr
%  Usage
%    ans = randmultivariatenormalabovethreshold(m,Sigma,thr);
%    ans = randmultivariatenormalabovethreshold(m,n,thr);
%  Inputs
%    m      number of replicates (number of vectors)
%    Sigma  n x n square covariance matrix; default is 1
%    n      size of vector, Sigma will be eye(n)
%    thr    threshold value; default is 0
%    maxit  number of iterations (optional)
%  Outputs
%    ans    random vector of size m times n
%  Description
%    This routine calls randnoutsideinterval
%  Note
%    This routine might need further verification
%  Examples
%  See also
%    help randnoutsideinterval


if (nargin == 0)
  m = 1;
  Sigma = 1;
  thr = 0;
elseif (nargin == 1)
  if (length(m) > 1)
     Sigma = m;
     m = 1;
     thr = 0;
  else
     Sigma = 1;
     thr = 0;
  end
elseif (nargin == 2)
  if (length(m) > 1)
     thr = Sigma;
     Sigma = m;
     m = 1;
  else
     thr = 0;
  end
end
if (length(Sigma) == 1) & Sigma(1,1) == round(Sigma(1,1)),
   n = Sigma;
   Sigma = eye(n);
else
   [n n2] = size(Sigma);
    if n~=n2,
       error('Input covariance matrix must be square')
    end
    if max(max(abs(Sigma-Sigma')))>eps,
       error('Input covariance matrix must be symmetric')
    end
    [E Lambda] = eig(Sigma);
    toll = 1.e-10;
    if any(diag(Lambda)< -toll),
       error('Input covariance matrix must be semi-positive definite')
    end
end
if nargin < 4,
   maxit = 10;
end

X = zeros(n,m);
stdev1 = sqrt(Sigma(1,1));
X(1,1:m) = stdev1*randnoutsideinterval(1,m,-thr/stdev1,thr/stdev1);
for i = 2:n,
   Sigma11 = Sigma(1:i-1,1:i-1);
   Sigma12 = Sigma(1:i-1,i);
   Sigma21 = Sigma(i,1:i-1);
   Sigma22 = Sigma(i,i);
   mui = Sigma21*(Sigma11\X(1:i-1,1:m));
   stdevi = sqrt(max(0,Sigma22 - Sigma21*(Sigma11\Sigma12)));
   X(i,1:m) = mui+stdevi*...
              randnoutsideinterval(1,m,(-thr-mui)/stdevi,(thr-mui)/stdevi);
end
if n > 1, for k = 1:maxit,
   for i = 1:n,
      noti = [(1:i-1) (i+1:n)];
      Sigma11 = Sigma(noti,noti);
      Sigma12 = Sigma(noti,i);
      Sigma21 = Sigma(i,noti);
      Sigma22 = Sigma(i,i);
      mui = Sigma21*(Sigma11\X(noti,1:m));
      stdevi = sqrt(max(0,Sigma22 - Sigma21*(Sigma11\Sigma12)));
      X(i,1:m) = mui+stdevi*...
          randnoutsideinterval(1,m,(-thr-mui)/stdevi,(thr-mui)/stdevi);
   end
end, end
   
ans = X;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
