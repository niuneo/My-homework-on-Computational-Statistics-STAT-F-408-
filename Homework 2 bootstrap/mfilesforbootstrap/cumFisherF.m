function p = cumFisherF(x,dof1,dof2);

% cumFisherF -- cumulative distribution function (cdf) for chi-square
%   Usage
%     p = cumFisherF(x,dof1,dof2);
%   Inputs
%     x     value in which cdf is evaluated
%     dof1  degrees of freedom (default 1)
%     dof2  degrees of freedom (default 1)
%   Outputs
%     p     cdf
%   Description
%     finds cumulative probability P(X < x) where X ~ F with (dof1,dof2) 
%     degrees of freedom.
%   See also
%     help cumgauss 
%     help cumchisquare
%   Note
%     cumFisherF is an alternative for FCDF if that function is not available

if nargin < 3,
   dof2 = 1,
   if nargin < 2,
      dof1 = 1;
   end
end

p = betainc(dof1*x./(dof1*x+dof2),dof1/2,dof2/2);
% Note that betainc is the NORMALISED incomplete beta function
