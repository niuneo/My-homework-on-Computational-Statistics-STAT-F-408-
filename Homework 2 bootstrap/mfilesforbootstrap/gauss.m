function gauss = gauss(x,s)

% gauss -- probability density function (pdf) for zero-mean normal RV
%   Usage
%     f = gauss(x,s);
%   Inputs
%     x     value in which pdf is evaluated
%     s     st.dev. of normal distribution (default 1)
%   Outputs
%     f     pdf
%   See also
%     help cumgauss (cumulative distribution function)
%     help Rcumgauss (survival function)
%     help multivargauss (multivariate normal density)


if (nargin == 1)
   s = 1;
end

gauss = exp(-x.*x./(2*s.*s))./(sqrt(2*pi)*s);

if length(s) == 1,
   if s<eps,
      gauss = zeros(size(x));
      gauss(abs(x)<eps) = Inf;
   end
else
   s0 = (s<eps);
   if any(s0),
      gauss(find(s0)) = 0;
      if length(x) > 1,
         x0 = (abs(x)<eps);
         s0x0 = s0.*x0;
         if any(s0x0),
            gauss(find(s0x0)) = Inf;
         end
      elseif abs(x)<eps,
         gauss(find(s0)) = Inf;
      end
   end
end


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
