
function T = exptimes(lambda,realrandom);

% usage: T = exptimes(lambda,realrandom);
% generates exponentially distributed data (time intervals) with parameter
% vector lambda
%
% See also: poisson.m (for poisson count generator) and randexp (for vector of
% exponentially distributed data with fixed lambda)

if nargin < 2 | ~realrandom,
   rand('seed',931316785);
end
lambda0 = (lambda < eps);
lambda = lambda + double(lambda0);
r = rand(size(lambda));
T = -log(1-r)./lambda;
T(find(lambda0)) = Inf;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
