
if exist('samplesize') ~= 1, samplesize = 100000; end
if exist('lambda') ~= 1, lambda = 1.4; end
if exist('alfa') ~= 1, alfa = 5.3; end
m = samplesize;
%n = samplesize;
X1 = randgammaMCMC(m,1,lambda,alfa);
X2 = randgamma(m,1,lambda,alfa);
% X1 = randgammaMCMC(n,1,lambda,alfa);
% X2 = randgamma(n,1,lambda,alfa);
plot(sort(X1),sort(X2),'.','color','b','linewidth',3)
axy = axis; b = max([X1;X2]);
hold on
plot([0 b],[0 b],'color','k','linewidth',3)
%hold off
axis([0 b 0 b])

%legend('randgammaMCMC','randgamma');
%set(gca,'fontsize',14,'fontweight','b')
xlabel('RVs by Metropolis-Hastings sampler')
ylabel('RVs by rejection sampling')
print -dpdf '3.pdf'
hold off;