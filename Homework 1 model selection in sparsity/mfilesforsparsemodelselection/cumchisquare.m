
function p = cumchisquare(x,dof);

% cumchisquare -- cumulative distribution function (cdf) for chi-square
%   Usage
%     p = cumchisquare(x,dof);
%   Inputs
%     x     value in which cdf is evaluated
%     dof   degrees of freedom (default 1)
%   Outputs
%     p     cdf
%   Description
%     finds cumulative probability P(X < x) where X ~ chi^2 with dof degrees of
%     freedom.
%   See also
%     help cumgauss 

if nargin < 2,
   dof = 1;
end
if dof < 1000,
   p = gammainc(x/2,dof/2);
   % note that Matlab's definition of the in-built function gammainc is
   % different from the standard definition: order of arguments is switched and
   % the function is normalised by the gamma of the second argument.
else
   p = cumgauss(x-dof,sqrt(2*dof));
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
