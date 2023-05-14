if exist('typeofX') ~= 1, typeofX = 'random'; end
if exist('realrandom') ~= 1, realrandom = 0; end
if exist('dimX') ~= 1, dimX = 0; end
   if ~dimX, m = 1000; n = 200; end % n = number of observations
clear dimX
if exist('snratio') ~= 1, snratio = 10; end
if exist('degreeofsparsity') ~= 1, degreeofsparsity = 0.05; end
if exist('bandwidth') ~= 1, bandwidth = 1; end
if exist('singvalues') ~= 1, singvalues = 1; end
if exist('n_singval') ~= 1, n_singval = 1; end
if exist('standardize') ~= 1, standardize = 0; end


p = degreeofsparsity;

[Y X beta stdev Z] = ...
  setupsimulationYisXbetaplussigmaZ(typeofX,m,n,p,snratio,realrandom,...
                                'bandwidth',bandwidth,...
                                'singularvalues',singvalues,...
                                'n_singval',n_singval);
mu = X*beta;
[n m] = size(X);

% iterative soft-thresholding with universal threshold
initit = 100;
maxit = 1000;
univthr = sqrt(2*log(n))*stdev;
betahats0 = iterativeST(Y,X,univthr,initit);
beta0 = betahats0;
MSEiterITSTuniv = zeros(1,maxit);
PEiterITSTuniv = zeros(1,maxit);
warning off
for it = 1:maxit,
   beta1 = iterativeST(Y,X,univthr,1,beta0);
   MSEiterITSTuniv(it) = norm(beta1-beta)^2/m;
   PEiterITSTuniv(it) = norm(X*beta1-mu)^2/n;
   beta0 = beta1;
end

% iterative ST for several thresholds
nthr = 50;
tt = (1:nthr)/nthr*univthr/5; % division by factor 5 after having observed that
                              % min PE error is way below univthr.
st = (1:nthr);

PEITST = zeros(1,nthr);
SSeITST = zeros(1,nthr);
k = zeros(1,nthr);
for t = 1:nthr,
   betahatITST = iterativeST(Y,X,tt(t),maxit,beta0);
    PEITST(t) = norm(X*betahatITST-mu)^2/n;

   k(t) = sum(abs(betahatITST)>eps);
   SSeITST(t) = norm(X*betahatITST-Y)^2/n;
%    PEtoSSe(t)=PEITST(t)./SSeITST(t)
   
end
CpITST = SSeITST+2*k*stdev^2/n-stdev^2;

% plot(st,PEtoSSe,'linewidth',3,'color','r')
% print -dpdf 'PE to SSE.pdf'

% plot(tt,PEITST,'linewidth',3,'color','r')
% 
% %a plot of the exact prediction error as a function of the number of iteration steps
% plot(st,PEITST,'linewidth',3,'color','r')
% title('Prediction error vs. number of iteration')
% xlabel('No. iteration')
% ylabel('PE')
% print -dpdf 'PE vs No.iteration.pdf'





% plot(tt,PEITST,'linewidth',3,'color','r')
% hold on
% [minPE opt] = min(PEITST);
% opthr = tt(opt);
% plot(tt(opt),minPE,'.','MarkerSize',10)
% str = ['(' num2str(tt(opt)) ',' num2str(minPE) ')'];
% text(tt(opt),minPE,str);
% title('Prediction error vs. soft threshold')
% xlabel('Threshold')
% ylabel('PE')
% print -dpdf 'PE vs soft threshold.pdf'
% hold off

plot(tt,PEITST,'linewidth',3,'color','r')
hold on
[minPE opt] = min(PEITST);
opthr = tt(opt);
plot(tt(opt),minPE,'.','MarkerSize',10)
str = ['(' num2str(tt(opt)) ',' num2str(minPE) ')'];
text(tt(opt),minPE,str);
title('PE*n/SD^2 vs. soft threshold')
xlabel('Threshold')
ylabel('PE*n/SD^2')
legend('Sample size =1000')
print -dpdf 'n=200 PE vs soft threshold.pdf'
hold off


plot(tt,CpITST,'linewidth',2,'color','b')
hold on
[minPE opt] = min(CpITST);
opthr = tt(opt);
plot(tt(opt),minPE,'.','MarkerSize',10)
str = ['(' num2str(tt(opt)) ',' num2str(minPE) ')'];
text(tt(opt),minPE,str);
title('Prediction error vs soft threshold')
xlabel('Threshold')
ylabel('Mallows Cp')
legend('Sample size =1000')
print -dpdf 'n=1000 Mallows Cp vs threshold.pdf'
hold off




plot(tt,PEITST,'linewidth',3,'color','r')
% xlabel('Threshold(?(2 lnn ) ?)/5')
% ylabel('PE')
hold on
plot(tt,CpITST,'linewidth',2,'color','b')
title('Sample size =1000')
xlabel('Threshold')
ylabel('Prediction error')
legend('PE*n/SD^2','studentized Cp')
print -dpdf 'n=1000 PE vs lamda Cp.pdf'
hold off


plot(k,PEITST,'linewidth',3,'color','r')
hold on
plot(k,CpITST,'linewidth',2,'color','b')
xlabel('No. nonzeros')
ylabel('PE, Cp')
legend('PE','Cp')
print -dpdf 'Cp, PE vs No nonzeros.pdf'
hold off

plot(k,PEITST,'linewidth',3,'color','r')
hold on
plot(k,CpITST,'linewidth',2,'color','b')
hold off

% iterative soft-thresholding with min PE threshold
betahats0 = iterativeST(Y,X,opthr,initit);
beta0 = betahats0;
MSEiterITSTopt = zeros(1,maxit);
PEiterITSTopt = zeros(1,maxit);
warning off
for it = 1:maxit,
   beta1 = iterativeST(Y,X,opthr,1,beta0);
   MSEiterITSTopt(it) = norm(beta1-beta)^2/m;
   PEiterITSTopt(it) = norm(X*beta1-mu)^2/n;
   beta0 = beta1;
end

plot(PEiterITSTopt,'.-','color','b')
hold on
plot(PEiterITSTuniv,'.-','color','r')
hold off

warning on

% debiasing
r = find(abs(betahatITST)>eps);
Xr = X(1:n,r);
betahatdeb = zeros(m,1);
betahatdeb(r) = (Xr'*Xr)\(Xr'*Y);

% n1 = sum(abs(betahatITST)>eps);
%n1 = sum(abs(beta)>eps)*3;

nthr = 50;
tt = (1:nthr)/nthr*univthr/5; % division by factor 5 after having observed that
                              % min PE error is way below univthr.
st = (1:nthr);

PEITST = zeros(1,nthr);
SSeITST = zeros(1,nthr);
k = zeros(1,nthr);
for t = 1:nthr,
   betahatITST = iterativeST(Y,X,tt(t),maxit,beta0);
    PEITST(t) = norm(X*betahatITST-mu)^2/n;

   k(t) = sum(abs(betahatITST)>eps);
   SSeITST(t) = norm(X*betahatITST-Y)^2/n;
%    PEtoSSe(t)=PEITST(t)./SSeITST(t)
   
end


n1 = m;
tt2 = (1:n1);
k_lars = zeros(1,m);
for t = 1:m,
   betahatLARS = LARS(Y,X,stdev,tt2(t)); 
   %iterativeST(Y,X,tt(t),maxit,beta0);
   XbetahatLARS(t);
   k_lars(t) = sum(abs(betahatLARS)>eps);
  
%    PEtoSSe(t)=PEITST(t)./SSeITST(t)
   
end

[betahatLARS XbetahatLARS xLARS CpLARS] = LARS(Y,X,stdev,200);

plot(betahatLARS,'linewidth',3,'color','r')
plot(CpLARS,'linewidth',3,'color','r')
print -dpdf 'CpLARS vs No nonzeros.pdf'
hold off



% plot(k_lars,PELARS,'linewidth',3,'color','r')
% hold on
plot(k_lars,CpLARS,'linewidth',2,'color','b')
xlabel('No. nonzeros')
ylabel('CpLARS')
legend('CpLARS')
print -dpdf 'CpLARS vs No nonzeros.pdf'
hold off

plot(betahatLARS,'linewidth',3,'color','r')
