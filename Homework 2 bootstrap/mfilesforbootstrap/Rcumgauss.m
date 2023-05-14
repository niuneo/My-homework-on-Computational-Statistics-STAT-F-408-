function Rcumgauss = Rcumgauss(x,s);

% Rcumgauss -- survival function --- complementary cumulative distribution 
%                                    function for zero-mean normal RV
%   Usage
%     f = Rcumgauss(x,s);
%   Inputs
%     x     value in which 1-cdf is evaluated
%     s     st.dev. of normal distribution (default 1)
%   Outputs
%     f     1-cdf
%   See also
%     help gauss (probability density function)
%     help cumgauss (cumulative distribution function)


if (nargin == 1)
   s = 1;
end

Rcumgauss = 0.5 * erfc(x./(sqrt(2)*s));

if length(s) == 1,
   if s<eps,
      Rcumgauss(abs(x)<eps) = 0;
   end
else
   s0 = (s<eps);
   if any(s0),
      if length(x) > 1,
         x0 = (abs(x)<eps);
         s0x0 = s0.*x0;
         if any(s0x0),
            Rcumgauss(find(s0x0)) = 0;
         end
      elseif abs(x)<eps,
         Rcumgauss(find(s0)) = 0;
      end
   end
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
