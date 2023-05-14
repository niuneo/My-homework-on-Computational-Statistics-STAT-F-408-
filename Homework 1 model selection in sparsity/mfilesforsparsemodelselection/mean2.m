
function M = mean2(X);

% ThreshLab/Noise/mean2 -- Computes mean of 1-d or 2-d data
%                         (vectors, signals; matrices, images)
%  Usage
%    M = mean2(X);
%  Inputs
%    X      observations
%  Outputs
%    M      scalar: overall mean
%  Description
%  See also
%    help var2


M = mean(mean(X));

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
