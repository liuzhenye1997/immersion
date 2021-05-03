%《Numerical_Optimization》中的L-BFGS
function [s,y,rho,d]=LBFGS(g,iteration,m,s,y,rho)
q=g;
rho(iteration)=1/(y(iteration,:)*s(iteration,:)');
alpha=zeros(iteration,1);
for i=iteration:-1:max(iteration-m+1,1)
    alpha(i)=rho(i)*s(i,:)*q;
    q=q-alpha(i)*y(i,:)';
end
gamma=s(iteration,:)*y(iteration,:)'/(y(iteration,:)*y(iteration,:)');
r=gamma*q;
for i=max(iteration-m+1,1):iteration
    beta=rho(i)*y(i,:)*r;
    r=r+s(i,:)'*(alpha(i)-beta);
end
d=r;