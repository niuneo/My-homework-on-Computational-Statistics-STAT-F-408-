
function x = invcumgauss(p,stdev);

% invcumgauss -- quantile function (inverse cdf) for zero-mean normal RV
%   Usage
%     x = invcumgauss(p,s);
%   Inputs
%     p     value in which quantile is evaluated, must be in [0,1]
%     s     st.dev. of normal distribution (default 1)
%   Outputs
%     x     quantile value
%   See also
%     help cumgauss (cdf = cumulative distribution function)


if nargin < 2, stdev = 1; end

if max(max(p))>1 | min(min(p))<0,
   error('invcumgauss: input must be probability values, i.e. between 0 and 1')
end
x = sqrt(2)*stdev*erfcinv(2*(1-p));

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
