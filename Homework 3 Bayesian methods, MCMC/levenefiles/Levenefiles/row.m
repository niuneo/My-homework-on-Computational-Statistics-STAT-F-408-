
function y = row(x);

% ThreshLab/General/row: makes sure that a vector has 1 row
%  Usage
%    y = row(x);
%  Inputs
%    x : row or column vector
%  Outputs
%    y : row vector
%  Description
%  Note
%  See also


if (size(x,1) > size(x,2) )
 y = x';
else
 y = x;
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
