function f = multivargauss(x,mu,Sigma)

% multivargauss -- probability density funct. (pdf) for multivariate normal RV
%   Usage
%     f = multivargauss(x,Sigma);
%     f = multivargauss(x,mu,Sigma);
%     f = multivargauss(x);
%   Inputs
%     x     vector or matrix (=row of vectors) in which pdf is evaluated
%     mu    vector or matrix (=row of vectors) of mean values
%     Sigma covariance matrix (must be the same for all vectors in x)
%   Outputs
%     f     pdf
%   Description
%     If x is a matrix, then each column of x is considered as a vactor in
%     which the evaluation takes place.
%     The mean vector can be different for every column in x, but the
%     covariance matrix must be the same for all columns in x
%   See also
%     help gauss (univariate, zero-mean normal density)

[n, m] = size(x);
if n == 1 & m > 1, x = x'; [n, m] = size(x); end
if nargin < 2,
   mu = zeros(size(x)); Sigma = eye(n); nmu = n; mmu = m;
else
   [nmu mmu] = size(mu);
   if nargin < 3 & mmu ~= m & mmu == nmu,
      Sigma = mu; mu = zeros(size(x));
   end
end

if nmu == 1 & mmu > 1, mu = mu'; [nmu mmu] = size(mu); end
if mmu == 1 & m > 1, mu = repmat(mu,1,m); mmu = m; end

x = x-mu;

f = exp(-sum(x.*(Sigma\x)/2))/sqrt(det(Sigma))/(2*pi)^(n/2);
