
function x = invcumabsgauss(p,mu,stdev,stopcrit);

% invcumabsgauss -- quantile function (inverse cdf) for absolute value of 
%                   normal RV
%   Usage
%     x = invcumabsgauss(p,mu,stdev);
%   Inputs
%     p     value in which quantile is evaluated, must be in [0,1]
%     mu    expected value of normal distribution (default 0)
%     stdev st.dev. of normal distribution (default 1)
%   Outputs
%     x     quantile value
%   See also
%     help invcumgauss


if nargin < 4, stopcrit = 1.e-6; end
if nargin < 3, stdev = 1; end
if nargin < 2, mu = 0; end

if max(max(p))>1 | min(min(p))<0,
   error('invcumgauss: input must be probability values, i.e. between 0 and 1')
end
x0 = Inf*ones(size(p));
x1 = invcumgauss(p/2+1/2,stdev)+abs(mu);
while max(max(abs(x1-x0)))/stdev > stopcrit,
   x0 = x1;
   % Apply Newton-Raphson iteration
   x1 = x0 - (cumgauss(x0-mu,stdev)-cumgauss(-x0-mu,stdev)-p)./...
             (gauss(x0-mu,stdev)+gauss(x0+mu,stdev));
end

a = isnan(x1);
x1(a) = invcumgauss(p(a),stdev)+abs(mu);
x = x1;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
