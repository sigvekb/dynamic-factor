function [xitt,xittm,Ptt,Pttm,loglik, Kgain]=KalmanFilter(initx,initV,x,A,C,R,Q)
% INPUTS
% x(:,t) - the observation at time t
% A - the system matrix
% C - the observation matrix 
% Q - the system covariance 
% R - the observation covariance
% initx - the initial state (column) vector 
% initV - the initial state covariance 
% OUTPUT:
% xittm = E[X(:,t) | y(:,1:t-1)]
% x values after transition update
% Pttm = Cov[X(:,t) | y(:,1:t-1)]
% Covariance matrix after transition update
% xitt = E[X(:,t) | y(:,1:t)]
% x values after observation update
% Ptt = Cov[X(:,t) | y(:,1:t)]
% Covariance matrix after observation update
% loglik - value of the loglikelihood

[T,n]=size(x);
rlag=size(A,1);

y=x';

% Initialization
xittm=[initx zeros(rlag,T)];
xitt=zeros(rlag,T);

Pttm=zeros(rlag,rlag,T);
Pttm(:,:,1)=initV;
Ptt=zeros(rlag,rlag,T);

logl = zeros(T,1);
Kgain = zeros(rlag,n,T);
% Forward pass over observed data
for j=1:T
    
    % See www.bzarg.com/p/how-a-kalman-filter-works-in-pictures/ for the
    % equations below
    L = C * Pttm(:,:,j) * C' + R;
    L_step = inv(chol(L)); 
    L_inv = L_step * L_step';
    K = (Pttm(:,:,j) * C') * L_inv; % / L;
    Kgain(:,:,j) = K;
    innovation = (y(:,j)-C*xittm(:,j));
    
    % Set K=0 for missing values
    [row, ~] = find(isnan(innovation));
    innovation(isnan(innovation)) = 0;
    K(:,row) = 0;
    
    % Update predictions after observation
    xitt(:,j) = xittm(:,j) + K * innovation;
    %Ptt(:,:,j)=Pttm(:,:,j) - K * C *Pttm(:,:,j);
    Ir = eye(rlag);
    Ptt(:,:,j) = (Ir - K*C) * Pttm(:,:,j) * (Ir - K*C)' + K * R * K'; % Based on Max Welling explanations
    
    % Get next transition predictions, predicting one-step-ahead
    xittm(:,j+1)= A * xitt(:,j);
    Pttm(:,:,j+1)= A * Ptt(:,:,j) * A' + Q;
    
    % Likelihood calculation not used
    % lik(j)=((2*pi)^(-N/2))*(abs((det(C*Pttm(:,:,j)*C'+R)))^(-.5))*...
    %    exp(-1/2*(y(:,j)-C*xittm(:,j))'*L*(-1/2*(y(:,j)-C*xittm(:,j))));
    
    e = y(:,j) - C*xittm(:,j); % error (innovation)
    e(isnan(e)) = 0;
    ss = length(A);
    d = size(e,1);
    S = C*Pttm(:,:,j)*C' + R;
    GG = C'*diag(1./diag(R))*C;
    
    % Mahalanobis distance calculations, find log likelihood
    detS = prod(diag(R))*det(eye(ss)+Pttm(:,:,j)*GG);
    denom = (2*pi)^(d/2)*sqrt(abs(detS));
    mahal = sum((e'/S) * e, 2);
    logl(j) = -0.5*mahal - log(denom);
end
loglik=sum(logl);
