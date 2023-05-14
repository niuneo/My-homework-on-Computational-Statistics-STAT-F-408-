
n = 2000;
m = 100;

a = -3;
b = 2;
stdev = 4;
mu = 5;

warning off
% you have to implement randCauchyplusNormal!
X = randCauchyplusNormal(m,n,mu,stdev,a,b);
Z = randn(m,n)*stdev+mu;
T = randcauchy(m,n)*b+a;
warning on

EXhat1 = mean(X,2);

IZ = 1./(1+(Z-a).^2/b^2);
IT = exp(-(T-mu).^2/(2*stdev^2));

EXhat2 = mean(Z.*IZ,2)./mean(IZ,2);
EXhat3 = mean(T.*IT,2)./mean(IT,2);


