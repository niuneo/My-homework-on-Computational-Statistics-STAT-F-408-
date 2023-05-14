nonzeros=zeros(1,50);
Cp=zeros(1,50);
for i=1:50,
    betal=iterativeST(Y,X,veclambda(i),500,beta0);
    iterativeST(Y,K,thr,maxit,betahat0)
    
    
    prederrorITST(i)=norm(X*betal-mu)^2/n;
    nnonzeros(i)=sum(abs(betal)>10*eps);
    dof=nnonzeros(i);
    muhat=X*betal;
    Cp(i)=norm(Y-muhat)^2/n+2*dof/n*stdev*2-stdev^2;
    1
end
    
   