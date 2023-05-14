
function xt = ST(x,t);

% ThreshLab/Threshold/ST -- soft-threshold function
%  Usage
%    xt = ST(x,t)
%  Inputs
%    x      input coefficients
%    t      threshold (scalar or vector of size(x))
%  Outputs
%    xt     soft-thresholded values of x
%  Description
%    ST(x,t) = (abs(x)-t).*sign(x).*(abs(x)>t)
%            = x.*(abs(x)>t).*(1-t./abs(x));
%  Note
%  See also:
%    help HT


xt = (x > t) .* (x - t) + (x < -t) .* (x + t);


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
