
function s2 = var2(X);

% ThreshLab/Noise/var2 -- variance of one or two-dimensional data
%                               (vectors, signals, matrices, images)
%  Usage
%    s2 = var2(X);
%  Inputs
%    X      observations
%  Outputs
%    s2     variance
%  Description
%    computes mean2(X.*X) - mean2(X)*mean2(X) = mean2((X-mean2(X)).^2)
%    no bias correction s2 = sum(X-mean(X))/n, not divided by n-1
%  See also
%    help mean2


s2 = mean2(X.*X) - mean2(X)*mean2(X);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
