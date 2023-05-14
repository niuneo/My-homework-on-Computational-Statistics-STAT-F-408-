function ratio = psnr(s,n)

% usage: ratio = psnr(s,n)
% computes psnr for signal s + n , n is noise and s is uncorrupted signal
ratio = 10 * log (max(abs(s))^2/var(n,1)) / log(10.0);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
