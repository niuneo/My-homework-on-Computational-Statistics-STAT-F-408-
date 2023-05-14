
n = 2000;
m = 100; % repeat 100 times estimations

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


figure; plot(EXhat1,'x'); hold on; 

plot(EXhat2,'or'); hold on;

plot(EXhat3,'c*');  hold on;

legend('theta1=EXhat1','theta2=EXhat2','theta3=EXhat3');
%set(gca,'fontsize',14,'fontweight','b')
xlabel('repeating times')
ylabel('expected value')
print -dpdf '1.pdf'
hold off;

var1=var(EXhat1,1)
var1=var(EXhat1,0)
var2=var(EXhat2,1)
var2=var(EXhat2,0)
var3=var(EXhat3,1)
var3=var(EXhat3,0)

xx=zeros(100,2);
xx(:,1)=EXhat1;
xx(:,2)=ones(100,1);

yy=zeros(100,2);
yy(:,1)=EXhat2;
yy(:,2)=2*ones(100,1);

zz=zeros(100,2);
zz(:,1)=EXhat3;
zz(:,2)=3*ones(100,1);

X1=[xx;yy;zz];


alpha=0.05;
Ltest=Levenetest(X1,alpha)