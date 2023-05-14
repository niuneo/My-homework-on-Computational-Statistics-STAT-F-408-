
function d = MAD1(y);

% MAD1 - ThreshLab/Noise -- median absolute deviation
%    Usage
%      d = MAD1(y);
%    Inputs
%      y: zero mean OR zero median vector (depending on how MAD is defined)
%    Outputs
%      d = median1(abs(y));
%    Description
%      Apply d = MAD1(y-median(y));
%        or  d = MAD1(y-mean(y));
%      depending on whether you want median absolute deviation from mean or 
%      median
%      Use
%        s = MAD1(y)/invcumgauss(0.75); 
%      to estimate noise deviation for zero-me(di)an sparse data.
%    See also
%      help empiricalquantile

d = median1(abs(y));

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
