
function [Y, X, beta, stdev, Z] = ...
    setupsimulationYisXbetaplussigmaZ(typeofdesign,m,varargin);

% setupsimulationYisXbetaplussigmaZ -- set up test with Y = X*beta + sigma*Z
%  Usage
%    [Y X beta stdev Z] =  setupsimulationYisXbetaplussigmaZ(typeofdesign,...
%                                               m,n,p,snratio,realrandom,...)
%  Inputs
%    typeofdesign possible models for design matrix. Currently available:
%             'eye' identity matrix (signal-plus-noise model: Y = beta+stdev*Z)
%             'random' random matrix
%             'randomtoeplitz' random selection of n rows from m x m 
%                              stationary filter matrix
%             'randombandtoeplitz' random selection of n rows from m x m 
%                                  bandlimited stationary filter matrix
%             'randomfromFourier' random selection of n rows from m x m Fourier
%                                 transformation matrix
%             'SVDgivenSV' U'*Sigma*V, with U and V random orthogonal matrices
%                          and Sigma = diag(sigma), where sigma is given as
%                          optional parameter in 'singularvalues'
%             'SVDrandomSV' U'*Sigma*V, with U and V random orthogonal matrices
%                          and Sigma = diag(sigma), where sigma r random
%                          singular values, where r given as optional parameter
%                          in 'n_singval'
%             'Vandermonde' Vandermonde matrix for given or random covariate
%                           values (declared as ,'vdmcovar',vdmcovar)
%    m      model size of beta
%    n      number of observations; default n = m
%    p      sparsity parameter in zero-inflated Laplacian model of beta:
%           p = P(beta =/= 0)      default p = 0.05
%    snratio signal-to-noise ratio; SNR = 20 log(norm(X*beta)/norm(noise))
%            default = 10;
%    realrandom  if false, then seed for random number generators is fixed
%                if true, then seed is current value of seed
%                default is false
%    'stdev',stdev use value stdev as standard deviation of nosie instead of
%                  signal-to-noise ratio
%    'normalizedcolumns',k normalize columns of X, using ell-k norm
%    'normalizedcolumns'   normalize columns of X, using ell-2 norm
%    'prior' or 'nonzerosmodel'   model for generating the nonzeros in the 
%                                 error-free data (beta)
%    examples:
%    'nonzerosmodel','Laplace'   zero-inflated Laplace for error-free data
%                                (default)
%    'prior','binary'    nonzeros are +1 or -1
%    'prior','constant'
%    
%  Outputs
%    Y = X*beta+stdev*Z
%    Y      observations, vector of size n x 1
%    X      design matrix, matrix of size n x m
%    beta   unobserved model, vector of size m x 1, generated according to zero
%           inflated Laplacian model with a = 1/5 and p from input
%    stdev  standard deviation of noise
%    Z      standardized observational noise, IID, normal vector of size n x 1
%  Examples of usage
%
%    use
%       beta = setupsimulationYisXbetaplussigmaZ('eye',n,n,p,Inf,realrandom);
%    to generate sparse vector only (no X, no errors)
%
%    if exist('typeofdesign') ~= 1, typeofdesign = 'random'; end
%    if exist('realrandom') ~= 1, realrandom = false; end
%    if exist('dimdesign') ~= 1, dimdesign = 0; end
%       if ~dimdesign, m = 1000; n = 200; end % n = number of observations
%    clear dimdesign
%    if exist('snratio') ~= 1, snratio = 10; end
%    if exist('degreeofsparsity') ~= 1, degreeofsparsity = 0.05; end
%    if exist('bandwidth') ~= 1, bandwidth = 1; end
%    if exist('singvalues') ~= 1, singvalues = min(m,n); end
%    if exist('n_singval') ~= 1, n_singval = 1; end
%    if exist('standardize') ~= 1, standardize = 0; end
%    if exist('vdmcovar') ~= 1, vdmcovar = sort(rand(1,n)); end
%
%    p = degreeofsparsity;
%
%    [Y,X,beta,stdev,Z] = setupsimulationYisXbetaplussigmaZ(typeofdesign,...
%                                               m,n,p,snratio,realrandom,...
%                                                  'bandwidth',bandwidth,...
%                                            'singularvalues',singvalues,...
%                                                  'n_singval',n_singval,...
%                                                    'vdmcovar',vdmcovar);
%    
% 

if nargin < 2,
   error('typeofdesign and m are mandatory inputs')
end
nnargin = nargin;
ar = 3; 
while ar < nargin, 
   if ischar(varargin{ar-2}), 
      nnargin = ar-1; 
      ar = nargin;
   end
   ar = ar+1;
end
if nnargin < 3, n = m; else n = varargin{1}; end
if nnargin < 4, p = 0.05; else p = varargin{2}; end
if nnargin < 5, snratio = 10; else snratio = varargin{3}; end
if nnargin < 6, realrandom = false; else realrandom = varargin{4}; end

nnargin = nnargin+1;
stdev0 = NaN;
vdmcovar = sort(rand(1,n));
r = min(m,n);
kerneltype = 'cos';
normalizedcolumns = NaN;
prior = 'Laplace';
for ar = nnargin:nargin, arg = varargin{ar-2}; if ischar(arg)
  switch arg
  case {'bandwidth'}
    bandwidth = varargin{ar-1};
  case {'singularvalues'}
    sigma = varargin{ar-1};
  case {'n_singval'}
    r = varargin{ar-1};
  case {'stdev'}
    stdev0 = varargin{ar-1};
  case {'vdmcovar'}
    vdmcovar = varargin{ar-1};
  case {'kernel','kerneltype','typekernel'}
    kerneltype = varargin{ar-1};
  case {'normalizedcolumns'}
    ell = varargin{ar-1};
    if isnumeric(ell),
       normalizedcolumns = ell;
    else
       normalizedcolumns = 2;
    end
  case {'prior','nonzerosmodel'}
    prior = varargin{ar-1};
  end
end, end
if typeofdesign(1:3) == 'SVD',
   if typeofdesign(4:6) == 'giv',
      r = length(sigma);
   else
      sigma = abs(rand(1,r));
   end
end
if length(typeofdesign)>6, if typeofdesign(2:7) == 'anderm',
   n = length(vdmcovar);
end, end

disp('-------------------------------------------------')
disp('Generating Y = X*beta+stdev*Z with following values:')
disp(['   typeofdesign = ' typeofdesign])
if ~isnan('normalizedcolumns'), if normalizedcolumns~= 0,
   disp(['             (columns of design (X) normalized, using ell-'...
                        num2str(normalizedcolumns) '-norm)'])
end, end
disp(['         n = ' num2str(n)])
disp(['         m = ' num2str(m)])
disp(['         p = ' num2str(p)])
disp(['       snr = ' num2str(snratio)])
disp(['realrandom = ' num2str(realrandom)])
disp('-------------------------------------------------')
disp(' '), disp(' ')

if realrandom==false,
   randn('state',2000);
   rand('state',2000);
end

Xeye = false;
switch typeofdesign,
case {'randomtoeplitz'}
   h = randn(1,m);
   X = toeplitz(h);
   r = randperm(m);
   r = sort(r(1:n));
   X = X(r,1:m);
case {'randombandtoeplitz'}
   h = randn(1,m);
   h(bandwidth+1:m) = 0;
   X = toeplitz(h);
   r = randperm(m);
   r = sort(r(1:n));
   X = X(r,1:m);
case {'random'}
   X = randn(n,m);
case {'randomuniform'}
   X = 2*rand(n,m)-1;
case {'randomband','bandrandom'}
   dotranspose = 0;
   if n>m, dotranspose = 1; nm = n; n = m; m = nm; end
   X = sparse(diag(randn(1,m)));
   for d = 1:bandwidth, 
       X = X+sparse(diag(randn(1,m-d),d))+sparse(diag(randn(1,m-d),-d));
   end
   r = [0 floor((1:m)/m*n)]; r = find(r(2:m+1)~=r(1:m));
   X = X(r,1:m);
   if dotranspose == 1, X = X'; nm = n; n = m; m = nm; end
case {'eye'}
   n = m;
   X = eye(m);
   Xeye = true;
case {'randomfromFourier'}
   halfm = floor(m/2-1/2);
   restm = m-1-halfm*2;
   X = fft(eye(m));
   X = X(1:m,1:1+halfm+restm);
   X = [real(X(:,1:1+halfm+restm)) imag(X(:,2:1+halfm))];
   X = X';
   r = randperm(m);
   r = sort(r(1:n));
   X = X(r,1:m);
case {'SVDgivenSV','SVDrandomSV'}
   mn = min(m,n);
   r = min(mn,r); sigma = sigma(1:r);
   sigma0 = zeros(1,mn); sigma0(1:r) = sigma;
   Sigma0 = diag(sigma0);
   Sigma = zeros(n,m); Sigma(1:mn,1:mn) = Sigma0;
   U = randn(n,n); V = randn(m,m);
   U = gramschmidt(U); V = gramschmidt(V);
   X = U'*Sigma*V;
case {'Vandermonde'}
   X = zeros(n,m);
   for k=1:m,
      X(1:n,k) = column(vdmcovar.^(k-1));
   end
case {'randomkernel'}
   center = (1:floor(m/n):m); center = center(1:n);
   for i = 1:n,
      o = center(i);
      beta = sort(rand(1,m));  betao = beta(o);
      X(i,1:m) = evalkernel(beta-betao,bandwidth,kerneltype);
      X(i,1:m) = X(i,1:m) /norm(X(i,1:m));
   end
end
if Xeye==false, % binary variable Xeye saves computation time in case X = eye(n)
   X = (1/norm(full(X)))*X;
end

if ~isnan(normalizedcolumns), if normalizedcolumns ~= 0,
   Xk = abs(X).^normalizedcolumns;
   columnnormsXk = (sum(Xk,1)).^(1/normalizedcolumns);
   X = X*diag(1./columnnormsXk);
end, end

a = 1/5;
r = rand(m,1);
beta = zeros(m,1);
nonzeros = find(r<p);
nnon0 = length(nonzeros);
switch prior
case {'Laplace','laplace'}
   beta(nonzeros) = randlaplace(nnon0,1,a);
case {'constant'}
   beta(nonzeros) = 1;
case {'binary'}
   beta(nonzeros) = round(rand(nnon0,1))*2-1;
end


if Xeye==false, Xbeta = X*beta; else, Xbeta = beta; end
if isinf(snratio), Y = Xbeta; else, Y = makenoisy(Xbeta,snratio); end
U = Y-Xbeta;
stdev = sqrt(var2(U));
if ~isnan(stdev0), U = U*stdev0/stdev; Y = Xbeta+U; end
Z = U/stdev;

disp('Simulation setup finished')
