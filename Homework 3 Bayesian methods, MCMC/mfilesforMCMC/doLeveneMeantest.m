
% see http://en.wikipedia.org/wiki/Levene's_test
%
% see doANOVA
% requires two input vectors of same size: Y and factor
% Y is the observational vector, factor contains for every observation a
% reference to the group to which the observation belongs. It should be an
% integer between 1 and k, with k the number of groups. If this is not the
% case, see doANOVA

Yhold = Y;
k = max(factor);
for p = 1:k,
   Yp = Y(find(factor==p));
   Ypmean(p) = mean(Yp);
end
residual = row(Y) - Ypmean(factor);
residual = reshape(residual,size(Y));
Y = abs(residual);
disp('Levene''s test statistic - Mean')
disp('------------------------------')
doANOVA
disp('------------------------------')
Z = Y;
Y = Yhold;
