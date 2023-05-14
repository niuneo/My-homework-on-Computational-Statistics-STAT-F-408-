
% requires vectors Y and factor; factor values should range from 1 to some k
% If this is not the case, then run the following code
%    factorvalues = unique(factor);
%    k = length(factorvalues);
%    for i = 1:n,
%       factor(i) = find(factorvalues==factor(i));
%    end
% Output variables: F0, Ypmean (mean of Y(factor==p)), np, SSWp, MSB, MSW,
% residual

[Levenetest] = Levenetest(X,alpha);

nY = length(Y);
k = max(factor);

% clear np Ypmean Yp SSWp
Ypmean = zeros(1,k);
np = zeros(1,k);
SSWp = zeros(1,k);
for p = 1:k,
   Yp = Y(find(factor==p));
   np(p) = length(Yp);
   Ypmean(p) = mean(Yp);
   SSWp(p) = sum((Yp-Ypmean(p)).^2);
end
residual = row(Y) - Ypmean(factor);
residual = reshape(residual,size(Y));
Ymean = mean(Y);
MSB = sum(np.*(Ypmean-Ymean).^2)/(k-1);
MSW = sum(SSWp)/(nY-k);
F0 = MSB/MSW
% F ~ F(k-1,nY-k) = F(4,304)
% F0 = 1.2580 is not significant
blomm = ((1:nY)-0.5)/nY;
plot(invcumgauss(blomm),sort(residual),'.')
% p = 1-FCDF(F0,k-1,nY-k)
p = 1-cumFisherF(F0,k-1,nY-k)
