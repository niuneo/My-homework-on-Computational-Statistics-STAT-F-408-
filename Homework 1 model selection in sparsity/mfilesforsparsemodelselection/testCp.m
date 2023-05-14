if exist('typeofK') ~= 1, typeofK = 'random'; end
if exist('realrandom') ~= 1, realrandom = 0; end
if exist('dimK') ~= 1, dimK = 0; end
   if ~dimK, m = 1000; n = 200; end % n = number of observations
clear dimK
if exist('snratio') ~= 1, snratio = 10; end
if exist('degreeofsparsity') ~= 1, degreeofsparsity = 0.05; end
if exist('bandwidth') ~= 1, bandwidth = 1; end
if exist('singvalues') ~= 1, singvalues = 1; end
if exist('n_singval') ~= 1, n_singval = 1; end
if exist('standardize') ~= 1, standardize = 0; end


p = degreeofsparsity;

% [Y X beta stdev Z] = ...
%   setupsimulationYisXbetaplussigmaZ(typeofX,m,n,p,snratio,realrandom,...
%                                 'bandwidth',bandwidth,...
%                                 'singularvalues',singvalues,...
%                                 'n_singval',n_singval);

[w K v stdev Z] = ...
  setupsimulationYisXbetaplussigmaZ(typeofK,m,n,p,snratio,realrandom,...
                                'bandwidth',bandwidth,...
                                'singularvalues',singvalues,...
                                'n_singval',n_singval);

% mu = X*beta;
% [n m] = size(X);

Kv = K*v;
[n m] = size(K);

% iterative soft-thresholding
initit = 100;
univthr = sqrt(2*log(n))*stdev;  % threshold 0~univthr
nthr = 50;
  thr = (1:nthr)/nthr*univthr/5;
 
prederrorITST = zeros(1,nthr);
CpITST = zeros(1,nthr);
for t = 1:nthr,
   vhats = iterativeST(w,K,thr(t),600);
   prederrorITST(t) = norm(K*vhats-Kv)^2/n;
   p1 = sum(abs(vhats)>eps)
   %CpITST(t) = norm(K*vhats-w)^2/n+2*p1/n*stdev^2-stdev^2;
CpITST(t) = (norm(K*vhats-w)^2+2*p1*stdev^2-n*stdev^2)/(n*stdev^2); % using the studentized version of Cp
end

% r=find(vhats==0)  % where betal==0, sparse
% plot(vhats,"r")

[minPE opt] = min(prederrorITST);
opthr = thr(opt);


plot(thr,prederrorITST,'r')
title('Prediction error vs. sample size')
xlabel('Thresholding')
ylabel('Prediction error')
hold on 
plot(thr,CpITST,'b')
hold off
legend('Studentized Cp')


plot(tt,PEITST,'linewidth',3,'color','r')
hold on
[minPE opt] = min(PEITST);
opthr = tt(opt);
plot(tt(opt),minPE,'.','MarkerSize',10)
str = ['(' num2str(tt(opt)) ',' num2str(minPE) ')'];
text(tt(opt),minPE,str);
title('Prediction error vs. soft threshold')
xlabel('Threshold')
ylabel('PE')
print -dpdf 'PE vs soft threshold.pdf'
hold off


% plot(thr,prederrorITST,'r')
% title('Prediction error vs. iteration')
% xlabel('Thresholding')
% ylabel('Prediction error')
% hold on 
% plot(thr,CpITST,'b')
% hold off
% legend('Studentized Cp')
