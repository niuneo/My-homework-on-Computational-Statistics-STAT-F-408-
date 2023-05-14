
alfa = 1/500;
beta = 1/50;
% alfa and beta are transition probabilities in a Markov Chain. 
% alfa = P(X(n) = 1| X(n-1) = -1) 
% beta = P(X(n) = -1| X(n-1) = 1)
% We want to create a MC where the state +1 lasts on average 1/alfa steps and
% the state -1 lasts on average 1/beta steps. 
% It can be proven that this corresponds to an Ising model with parameters tau
% and gamma (drift) equal to:
gamma = log((1-alfa)/(1-beta))/2
tau = log(alfa*beta/((1-alfa)*(1-beta)))/4
kappa = 0.3;
K = 4;
sigma = 1;

p = 10000;
X = randising(p,1,tau,gamma);
stdevM = kappa*(X==-1) + K*(X==1);
M = randn(size(X)).*stdevM;
V = randn(size(X))*sigma;
Y = M+V;

% Gibbs sampling

% initialisation

thruniv = sqrt(2*log(p))*sqrt(kappa^2+sigma^2);
Xuniv = (abs(Y)>thruniv)*2-1;

% We first sample the prior to find the prior marginal probabilities

Xn = randising(Xuniv,tau,gamma,10);
XX = Xn;
nsteps = 100;
for n = 1:nsteps-1,
    Xn = randising(Xn,tau,gamma,1);
    XX = XX+Xn;
end
XX = XX/nsteps;

