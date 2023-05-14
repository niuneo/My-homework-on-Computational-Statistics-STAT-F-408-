
function y = evalkernel(x,h,kerneltype)

% evalkernel - ThreshLab2/Irregular/MultiscaleKernel/ --- evaluation of kernel 
%              in kernel smoothing or kernel density estimation
%   Usage
%      y = evalkernel(x,h,kerneltype)
%   Inputs
%      x      x-coordinates in which kernel is evaluated
%      h      kernel bandwidth
%      kerneltype  string-type variable, alternatives are
%           'none','dirac','Dirac','kronecker','Kronecker'
%           'uniform','Uniform','unif'
%           'Triangle','triangle','tri'
%           'Epanechnikov'
%           'Quartic','quartic'
%           'Triweight','triweight'
%           'cosine','Cosine','cos'
%           'Gaussian','gaussian','gauss','Gauss'
%   Outputs
%      y    function values in x
%   See also
%      help kernelestimation


% x/h is supposed to be on [-1,1] or ]-Inf,Inf[
% kernels have integral equal to bandwidth h

if nargin < 2,
   h = 1;
end
if nargin < 3,
   if ischar(h),
      kerneltype = h; h = 1;
   else
      kerneltype = 'uniform';
   end
end

if length(h) == 1,
   x = x/h;
elseif length(h) == length(x),
   if size(h)~=size(x), h = h'; end
   x = x./h;
end
y = (abs(x) <= 1);

switch kerneltype
case {'none','dirac','Dirac','kronecker','Kronecker'}
   y = zeros(size(x))/2; y(x==0) = 1;
case {'uniform','Uniform','unif'}
   y = ones(size(x))/2 .* y;
case {'Triangle','triangle','tri'}
   y = (1 - abs(x)) .*y;
case {'Epanechnikov'}
   y = 3*(1-x.^2)/4.*y;
case {'Quartic','quartic'}
   y = 15*(1-x.^2).^2/16 .*y;
case {'Triweight','triweight'}
   y = 35*(1-x.^2).^3/32 .*y;
case {'cosine','Cosine','cos'},
   y = pi*cos(pi*x/2)/4 .* y;
case {'Gaussian','gaussian','gauss','Gauss'}
   y = gauss(x);
end
