function  ratio = snr2(s,n)

% ThreshLab/Noise/snr2 -- Computes signal-to-noise ratio for 2-dimensional 
%                         data (matrices, images)
%  Usage
%    ratio = snr2(s,n);
%  Inputs
%    s      Data without noise
%    n      noise
%  Outputs
%    ratio  signal-to-noise ratio
%  Description
%    computes snr for signal s + n , n is noise and s is uncorrupted data
%    ratio = 10 log10(var(s)/var(n))
%  See also
%    help snr

ratio = 10 * log (var2(s)/var2(n)) / log(10.0);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
