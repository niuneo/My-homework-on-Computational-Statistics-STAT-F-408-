
function y = makenoisy(f,B,realrandom,ratio);

% makenoisy -- adds normal noise to data
%  Usage
%    y = makenoisy(f,B,realrandom,ratio);
%    y = makenoisy(f,ratio);
%  Inputs
%    f          Data without noise
%    B          filter for colored/correlated noise (B=1 for white/independent
%               noise)
%    realrandom binary value: 0 (default) for noise with fixed seed
%                             1 for non-repeatable real noise
%    ratio      signal-to-noise ration = 10 log10(var(f)/var(noise))
%  Outputs
%    y          Data with noise
%  Description
%    set B = 1 for white noise.
%    for colored noise with given frequency spectrum, use signal processing
%    toolbox: B = fir1(100,0.6,'high'); for a high-pass filter (more info on
%    help fir1)
%  See also
%    help snr
%
% reported problem: not yet upgraded to work with non-trivial 2-D entry for B

if nargin == 2,
   ratio = B;
   realrandom = 0;
   B = 1;
end
if nargin == 3,
   if length(B) == 1,
      ratio = realrandom;
      realrandom = B;
      B = 1;
   else
      ratio = realrandom;
      realrandom = false;
   end
end
if (realrandom == false)
  randn('seed',931316785);
end
whitenoise = randn(size(f));
if (length(B) > 1)
   noise = filter(B,1,whitenoise);
else
   noise = whitenoise;
end

% Following lines to obtain desired noise level
snr1 = snr2(f,noise);
factor = 10^((snr1-ratio)/20);
noise = noise*factor;

y = f+noise;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
