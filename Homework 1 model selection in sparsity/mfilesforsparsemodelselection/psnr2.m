function  ratio = psnr2(s,n)
% usage: ratio = psnr2(s,n)
% computes psnr for image s + n , n is noise and s is uncorrupted image
ratio = 10 * log (max(max(abs(s)))^2/var2(n)) / log(10.0);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
