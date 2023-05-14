function ratio = snr(s,n)

% ThreshLab/Noise/snr -- Computes signal-to-noise ratio
%  Usage
%    ratio = snr(s,n);
%  Inputs
%    s      Data without noise
%    n      noise
%  Outputs
%    ratio  signal-to-noise ratio
%  Description
%    computes snr for signal s + n , n is noise and s is uncorrupted signal
%    ratio = 10 log10(var(s)/var(n))
%  See also
%    help snr2

ratio = 10 * log (var(s,1)/var(n,1)) / log(10.0);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
