function [xitT,PtT,PtTm]=KalmanSmoother(A,xitt,xittm,Ptt,Pttm,C,R)
% INPUTS
% y(:,t) - the observation at time t
% A - the system matrix
% xittm = E[X(:,t) | y(:,1:t-1)]
% Pttm = Cov[X(:,t) | y(:,1:t-1)]
% xitt = E[X(:,t) | y(:,1:t)]
% Ptt = Cov[X(:,t) | y(:,1:t)]
% C - the observation matrix 
% R - the observation covariance

% OUTPUT:
% xitT = E[X(:,t) | y(:,1:T)]
% PtT = Cov[X(:,t) | y(:,1:T)]
% PtTm = Cov[X(:,t+1),X(:,t) | y(:,1:T)]

[T]=size(xitt,2);
r=size(A,1);
Pttm=Pttm(:,:,1:end-1);
xittm=xittm(:,1:end-1);
J=zeros(r,r,T);


for i=1:T-1
    J(:,:,i)= (Ptt(:,:,i) * A') / Pttm(:,:,i+1);
end

[var, factors] = size(C);
L = zeros(var, var, T);
K = zeros(factors, var, T);
for i=1:T
    L(:,:,i)= C * Pttm(:,:,i) * C' + R;
    K(:,:,i)= (Pttm(:,:,i) * C') / L(:,:,i);
end


xitT=[zeros(r,T-1)  xitt(:,T)];
PtT=zeros(r,r,T);
PtTm=zeros(r,r,T);
PtT(:,:,T)=Ptt(:,:,T);
PtTm(:,:,T)=(eye(r)-K(:,:,T)*C)*A*Ptt(:,:,T-1);

for j =1:T-1
    
    xitT(:,T-j)= xitt(:,T-j)+J(:,:,T-j)*(xitT(:,T+1-j)-xittm(:,T+1-j));
    PtT(:,:,T-j)=Ptt(:,:,T-j)+J(:,:,T-j)*(PtT(:,:,T+1-j)-Pttm(:,:,T+1-j))*J(:,:,T-j)';
    
    
end


for j =1:T-2
    PtTm(:,:,T-j)=Ptt(:,:,T-j)*J(:,:,T-j-1)'+J(:,:,T-j)*(PtTm(:,:,T-j+1)-A*Ptt(:,:,T-j))*J(:,:,T-j-1)';
end