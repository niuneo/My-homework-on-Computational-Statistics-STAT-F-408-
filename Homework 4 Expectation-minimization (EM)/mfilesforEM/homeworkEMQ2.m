
if exist('observesign') ~= 1, observesign = false; end
n = 1000;
Z = randn(n,1);
mu = 2;
stdev = 3;
% stdev = 0.1;
X = mu+stdev*Z;
m = 100;
mu0 = 0;
mu1 = 3;
muvec = (0:m)/m*(mu1-mu0)+mu0;

if observesign == true, Y = X; else Y = abs(X); end

[LL stdev2hat] = profileloglikelihoodabsnormalfixedpoint(muvec,X);
stdevhat = sqrt(stdev2hat);

% now implement LLX = profile likelihood, X observed including sign(X), as
% function of muvec
% (+/- 5 lines)
%

figure(2)
subplot(2,1,1)
plot(muvec,LL,'r')
legend('LLprofile(\mu) for |X|','location','NW')
subplot(2,1,2)
plot(muvec,LLX)
hold off
legend('LLprofile(\mu) for X','location','NW')

smax = 5; smin = 0;
s = (0:m)/m*(smax-smin)+smin;
s2 = s.^2;
s = repmat(s,n,1);
YY = repmat(Y,1,m+1);
LLabsXmuknown = sum(log(gauss(YY-mu,s)+gauss(YY+mu,s)),1);
LLabsXmuis1 = sum(log(gauss(YY-1,s)+gauss(YY+1,s)),1);

figure(3)
plot(s2,LLabsXmuknown)
hold on
plot(s2,LLabsXmuis1,'r')
hold off
axis([0 25 -2120 -1980])
legend('LL(\sigma^2) for true \mu','LL(\sigma^2) for \mu = 1')

% interprete these 2 curves

[muhatEM stdev2hatEM itEM] = EMabsnormalRV(X);
