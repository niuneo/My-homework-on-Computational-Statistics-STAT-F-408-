
function y = column(x)

% usage: y = column(x);
% input is row vector or column vector x
% output is column vector y
% Maarten Jansen 7 november 1997

if (size(x,1) < size(x,2) )
 y = x';
else
 y = x;
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
