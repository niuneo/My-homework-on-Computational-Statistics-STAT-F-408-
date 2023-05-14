
function m = median1(X);

% usage: m = median1(X);
% computes global median of data X. X may be 1D or 2D.

m = median(reshape(X,size(X,1)*size(X,2),1));

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
