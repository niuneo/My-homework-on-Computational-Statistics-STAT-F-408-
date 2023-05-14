
% X ~ Gamma(lambdaX,alfaX) Y ~ Gamma(lambdaY,alfaY)
% theta = E(X)/E(Y)

% rand('state',0);

alfaX = 4.3;
alfaY = 5.7;
lambdaX = 1;
lambdaY = 0.8;
theta = (alfaX*lambdaY)/(alfaY*lambdaX);

Nsimul = 100000;
if exist('nsample') ~= 1, nsample = 10; end
nsample = 500;
B = 10000;

X = randgamma(1,nsample,lambdaX,alfaX);
Y = randgamma(1,nsample,lambdaY,alfaY);
R = X./Y;

%%%%% non-parametric bootstrap %%%%%%%

rand('state',0);

T = mean(X)/mean(Y);
r1 = ceil(rand(B,nsample)*nsample);
Xstar = X(r1);
r2 = ceil(rand(B,nsample)*nsample);
Ystar = Y(r2);
Tstar = mean(Xstar,2)./mean(Ystar,2);

%%%%% parametric bootstrap %%%%%%%
% We assume that gamma model of X and Y are known

Xbar = mean(X);
Ybar = mean(Y);
Xsquaredbar = mean(X.^2);
Ysquaredbar = mean(Y.^2);

% moment estimators of parameters
alfaXhat = Xbar^2/(Xsquaredbar-Xbar^2);
lambdaXhat = Xbar/(Xsquaredbar-Xbar^2);
alfaYhat = Ybar^2/(Ysquaredbar-Ybar^2);
lambdaYhat = Ybar/(Ysquaredbar-Ybar^2);

XstarP = randgamma(B,nsample,lambdaXhat,alfaXhat);
YstarP = randgamma(B,nsample,lambdaYhat,alfaYhat);
TstarP = mean(XstarP,2)./mean(YstarP,2);


%%%%%% We can simulate the true distribution of the statistic

Xsimul = randgamma(Nsimul,nsample,lambdaX,alfaX);
Ysimul = randgamma(Nsimul,nsample,lambdaY,alfaY);
Tsimul = mean(Xsimul,2)./mean(Ysimul,2);

%%%%%% We compare the density functions

M = 1000;
tmin = min([min(Tstar),min(TstarP),min(Tsimul)]);
tmax = max([max(Tstar),max(TstarP),max(Tsimul)]);
tt = tmin + (0:M)/M * (tmax-tmin);
h = (tt(M)-tt(1))/10;

fTstar  = kernelestimation(tt,Tstar,h,'cos');
fTstarP = kernelestimation(tt,TstarP,h,'cos');
fTsimul = kernelestimation(tt,Tsimul,h,'cos');

% plot(tt,fTstar,'g','linewidth',2)
% hold on
% plot(tt,fTstarP,'r--','linewidth',2)
% plot(tt,fTsimul,'b','linewidth',2)
% hold off
% legend('Non-param. bootstr.','Param. bootstrap','True distr. f_T(t)')
% set(gca,'fontsize',14,'fontweight','b')

%%%%%%%%%variance question c
biasT = (alfaX/alfaY)*(lambdaY/lambdaX)*(1./(1-1./(nsample*alfaY))-1);
biasT2 = mean(Tsimul)- theta; 

biasTstar= mean(Tstar)-T;
biasTstarP= mean(TstarP)-T;


Tstar = mean(Xstar,2)./mean(Ystar,2);
TstarP = mean(Xstar,2)./mean(Ystar,2);

varTsimul = var2(Tsimul);
varTstar = var2(Tstar);
varTstarP = var2(TstarP);

%%%%%%%%question d, confidence interval
tau = median(Tsimul); 
delta = tau - theta;
tau_starP = median(TstarP); 
%delta_starP = tau_starP - mean(TstarP);
delta_starP=tau_starP-(alfaXhat*lambdaYhat)/(alfaYhat*lambdaXhat);

tau_star = median(Tstar); 
% delta_star = tau_star - mean(Tstar);
delta_star = tau_star - T;

%%%%%% calculate CI
alfa = 0.05; 
Tstarsorted = sort(Tstar); 
TstarPsorted = sort(TstarP);
i = round([alfa/2 1-alfa/2]*B);
percentilebootstrapCInonpar = Tstarsorted(i)';
percentilebootstrapCIpar = TstarPsorted(i)';
correctedpercentilebootstrapCInonpar=percentilebootstrapCInonpar+delta_star;
correctedpercentilebootstrapCIpar=percentilebootstrapCIpar+delta_starP;
%%%%%% calculate basic CI, question e
% standard error SE(theta) by bootstrapping
% inds=unidrnd(nsample, nsample, B);
% thetaboot=R(inds);
% T2 = mean(thetaboot,1)
% SEtheta2=std(T2);
%SEtheta=std(Tstar)/sqrt(length(Tstar));
 
SEtheta=std(Tstar); % Compute Standard Error
CI95 = tinv([alfa/2 1-alfa/2], length(Tstar)-1);    % Calculate 1-alpha Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, SEtheta, CI95(:));           % Calculate 1-alpha Confidence Intervals 
CIbasic=yCI95+mean(Tstar);     % Calculate basic Intervals

%%%%%%%%question f
%Part I
Xbar = mean(X); 
Ybar = mean(Y);
T = Xbar/Ybar;
Xstarbar = mean(Xstar,2);  % Xstarbar is B means by bootstrapping; Xstar is a matrix of B lines and nsample columns
Ystarbar = mean(Ystar,2); % Ystarbar is B means by bootstrapping; Ystar is a matrix of B lines and nsample columns

SX2 = sum((X-Xbar).^2)/(nsample-1); % calculate S_X^2  
SY2 = sum((Y-Ybar).^2)/(nsample-1);% calculate S_Y^2 
SXbar2 = SX2/nsample; % calculate S_Xbar^2 (Calculate Standard Error Of The Mean)
SYbar2 = SY2/nsample; % calculate S_Ybar^2 (Calculate Standard Error Of The Mean)
A = SXbar2/(Xbar^2-SXbar2)+SYbar2/(Ybar^2-SYbar2);
ST2 = T^2*A/(A+1); % an estimator for S_T^2

%Part II, an advanced version. The bootstrap version does exactly the same, but with stars. 
SXstar2 = sum((Xstar-Xstarbar*ones(1,nsample)).^2,2)/(nsample-1); % resamples derived from B bootstrapping replicates
SYstar2 = sum((Ystar-Ystarbar*ones(1,nsample)).^2,2)/(nsample-1); % resamples derived from B bootstrapping replicates
SXstarbar2 = SXstar2/nsample; % Calculate Standard Error Of The Mean
SYstarbar2 = SXstar2/nsample; % Calculate Standard Error Of The Mean
Astar = SXstarbar2./(Xstarbar.^2-SXstarbar2)+...
        SYstarbar2./(Ystarbar.^2-SYstarbar2);
STstar2 = Tstar.^2.*Astar./(Astar+1); % an estimator for S_T*^2

%Part III
% SEstar(b) is the estimated standard error of Tstar(b)
           
%calculate "zvals= [Tstar(b) -T]/SEstar(b)", where Tstar(b) is the boostrap replicate of T 
zvals_simple=(Tstar-T)./sqrt(ST2);      
zvals_advanced=(Tstar-T)./sqrt(STstar2);  
alfa = 0.05; 
k=B*alfa/2;
%simple version
    %get the quantiles
    szval_simple = sort(zvals_simple); 
    tlo_simple=szval_simple(k);
    thi_simple=szval_simple(B-k);
    %get the endpoints of the interval
    blo_simple=T-thi_simple* sqrt(ST2);
    bhi_simple=T-tlo_simple*sqrt(ST2);
%advanced version
   %get the quantiles
   szval_advanced = sort(zvals_advanced); 
   tlo_advanced=szval_advanced(k);
   thi_advanced=szval_advanced(B-k);
   %get the endpoints of the interval
   blo_advanced=T-thi_advanced* sqrt(mean(STstar2));
   bhi_advanced=T-tlo_simple*sqrt(mean(STstar2));
   
alfa = 0.05; 
Tstarsorted = sort(Tstar); 
i = round([alfa/2 1-alfa/2]*B);
percentilebootstrapCInonpar = Tstarsorted(i)';
percentilebootstrapCIpar = TstarPsorted(i)';
correctedpercentilebootstrapCInonpar=percentilebootstrapCInonpar+delta_star;
correctedpercentilebootstrapCIpar=percentilebootstrapCIpar+delta_starP;

%%%%%%bias estimates
biasT = (alfaX/alfaY)*(lambdaY/lambdaX)*(1./(1-1./(nsample*alfaY))-1)
biasTstar = (alfaX/alfaY)*(lambdaY/lambdaX)*(1./(1-1./(nsample*alfaY))-1)
biasTstarP = (alfaX/alfaY)*(lambdaY/lambdaX)*(1./(1-1./(nsample*alfaY))-1)

%%%%%%% Regression problem

% x1 = [-2 -1 0 1 2]';
% x2 = [-3 -2 1 2 3]';
% n = length(x1);
% y = 3+5*x1+4*x2-2*x1.*x2.^2;
% alfa = 2; lambda = sqrt(2);
% epsilon = randgamma(size(x2),lambda,alfa) - alfa/lambda;
% X = [ones(size(x1)) x1 x2 x1.*x2.^2];
% beta = [3 5 4 -2]';
% y = X*beta;
% Y = y+epsilon;
% betahat = (X'*X)\(X'*Y)
