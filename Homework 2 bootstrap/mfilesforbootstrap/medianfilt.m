
function y = medianfilt(x,taps)

% medianfilt -- apply median filter
%  Usage
%    y = medianfilt(x,taps)
%  Inputs
%    x      input data
%    taps   number of taps (=neighbors in filter)
%  Outputs
%    y      filtered data
%  Description
%    finds running median of (x-taps/2:x+taps/2)
%  Note
%  See also

if nargin < 2,
   taps = 3;
end
transpose = 0;
if size(x,2) == 1,
   x = x';
   transpose = 1;
end
Y = row(x);
N = length(x);
M = taps;
Y = repmat(x,M,1);
m = floor(M/2);
for r=2:M
   Y(r,r:N) = x(1:N-r+1);
end
y = shift(median(Y),m);
if transpose==1, y = y'; end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
