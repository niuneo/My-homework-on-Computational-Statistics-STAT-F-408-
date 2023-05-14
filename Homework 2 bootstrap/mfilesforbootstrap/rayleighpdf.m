
function f = rayleighpdf(x,s);

% ThreshLab/Noise/rayleighpdf -- Rayleigh probablity density function 
%    Usage
%      f = rayleighpdf(x,s);
%    Inputs
%      x    vector or scalar
%      s    parameter vector (of size(x) if size(x)>1) or scalar
%    Outputs
%      f    vector of size(x) and/or size(s)
%    Description
%      f = x/s^2*exp(-x/(2*s^2))

f = x./s^2.*exp(-x./(2*s.^2));

if length(s) == 1,
   if s<eps,
      f = zeros(size(x));
      f(abs(x)<eps) = Inf;
   end
else
   s0 = (s<eps);
   if any(s0),
      f(find(s0)) = 0;
      if length(x) > 1,
         x0 = (abs(x)<eps);
         s0x0 = s0.*x0;
         if any(s0x0),
            f(find(s0x0)) = Inf;
         end
      elseif abs(x)<eps,
         f(find(s0)) = Inf;
      end
   end
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 

