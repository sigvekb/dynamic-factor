function [beta_AQ, delta_C, gamma1_AQ, gamma2_AQ, ...
          x1, V1, loglik_t, xsmooth, delta, gammaCmiss, gamma] = ...
    Estep(y, A, C, Q, R, initx, initV, block, W)

% This function computes the (expected) sufficient statistics for a single Kalman filter sequence.
%
% y is the observable and x the hidden state

% INPUTS
% y(:,t) - the observation at time t
% A - the system matrix
% C - the observation matrix 
% Q - the system covariance 
% R - the observation covariance
% initx - the initial state (column) vector 
% initV - the initial state covariance 

% OUTPUTS: the expected sufficient statistics, i.e.
% beta = sum_t=1^T (x_t * x'_t-1)
% gamma = sum_t=1^T (x_t * x'_t)  
% delta = sum_t=1^T (y_t * x'_t)
% gamma1 = sum_t=2^T (x_t-1 * x'_t-1)
% gamma2 = sum_t=2^T (x_t * x'_t)
% x1  expected value of the initial state
% V1  variance of the initial state
% loglik value of the loglikelihood
% xsmooth expected value of the state

[n, T] = size(y);

% Use the Kalman smoother to compute 
% xsmooth = E[X(:,t) | y(:,1:T)]
% Vsmooth = Cov[X(:,t) | y(:,1:T)]
% VVsmooth = Cov[X(:,t), X(:,t-1) | y(:,1:T)] t >= 2
% loglik = sum{t=1}^T log P(y(:,t))
[xitt,xittm,Ptt,Pttm,loglik_t]=KalmanFilter(initx,initV,y',A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=KalmanSmoother(A,xitt,xittm,Ptt,Pttm,C,R);

% Create all versions of xsmooth, y, Vsmooth, VVsmooth needed
% Global + block
[r, ~] = size(xsmooth);

% All factors
xSmooth_AQ = cell(1, r);

for i=1:r
    m = xsmooth(i,:);
    xSmooth_AQ{i} = m;
end

Vsmooth_AQ = cell(1, r);
for i=1:r
    V = Vsmooth(i,i,:);
    Vsmooth_AQ{i} = V;
end

VVsmooth_AQ = cell(1, r);
for i=1:r
    VV = VVsmooth(i,i,:);
    VVsmooth_AQ{i} = VV;
end

gamma_AQ = zeros(1, r);
beta_AQ = zeros(1, r);

for t=1:T
    % AQ
    for i=1:r
        gamma_AQ(i) = gamma_AQ(i) + xSmooth_AQ{i}(:,t)*xSmooth_AQ{i}(:,t)' + Vsmooth_AQ{i}(:,:,t);
        if t>1
            beta_AQ(i) = beta_AQ(i) + xSmooth_AQ{i}(:,t)*xSmooth_AQ{i}(:,t-1)' + VVsmooth_AQ{i}(:,:,t); 
        end
    end
    
end
gamma1_AQ = zeros(1, r);
gamma2_AQ = zeros(1, r);
for i=1:r
    gamma1_AQ(i) = gamma_AQ(i) - xSmooth_AQ{i}(:,T)*xSmooth_AQ{i}(:,T)' - Vsmooth_AQ{i}(:,:,T);
    gamma2_AQ(i) = gamma_AQ(i) - xSmooth_AQ{i}(:,1)*xSmooth_AQ{i}(:,1)' - Vsmooth_AQ{i}(:,:,1);
end

%-------------------------------------------------------------------------
% Values without consideration of block-structure
delta = zeros(n,r);
gamma = zeros(r,r);
gammaCmiss = zeros(n*r, n*r);
beta = zeros(r, r);
y_temp = y;
y_temp(isnan(y_temp)) = 0;
for t=1:T
    deltaChange = W(:,:,t)*y_temp(:,t)*xsmooth(:,t)';
    delta = delta + deltaChange; % D
    gammaChange = xsmooth(:,t)*xsmooth(:,t)' + Vsmooth(:,:,t);
    missingValueKron = kron(gammaChange, W(:,:,t));
    gammaCmiss = gammaCmiss + missingValueKron; % C
    gamma = gamma + gammaChange;
    if t>1
        beta = beta + xsmooth(:,t)*xsmooth(:,t-1)' + VVsmooth(:,:,t); % 
    end
end
gamma1 = gamma - xsmooth(:,T)*xsmooth(:,T)' - Vsmooth(:,:,T); % 
gamma2 = gamma - xsmooth(:,1)*xsmooth(:,1)' - Vsmooth(:,:,1);
%-------------------------------------------------------------------------

x1 = xsmooth(:,1);
V1 = Vsmooth(:,:,1);

