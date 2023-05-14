
function y = row(x);

% function y = row(x);
% input is column-vector or row-vector
% output is rowvector
% Maarten Jansen 14 november 1997

if (size(x,1) > size(x,2) )
 y = x';
else
 y = x;
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
