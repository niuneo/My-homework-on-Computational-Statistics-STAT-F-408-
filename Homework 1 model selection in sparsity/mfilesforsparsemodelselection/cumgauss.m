function cumgauss = cumgauss(x,s);

% cumgauss -- cumulative distribution function (cdf) for zero-mean normal RV
%   Usage
%     f = cumgauss(x,s);
%   Inputs
%     x     value in which cdf is evaluated
%     s     st.dev. of normal distribution (default 1)
%   Outputs
%     f     cdf
%   See also
%     help gauss (probability density function)
%     help Rcumgauss (survival function)
%     help invcumgauss (quantile function)

if (nargin == 1)
   s = 1;
end

% following expression numerically unstable for x -> -Inf
% cumgauss = 1-Rcumgauss(x,s);
cumgauss = Rcumgauss(-x,s);


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
