function [beta_AQ, gamma_C, delta_C, gamma1_AQ, gamma2_AQ, x1, V1, loglik_t, xsmooth] = ...
    Estep(y, A, C, Q, R, initx, initV, block)

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

[~, T] = size(y);

% Use the Kalman smoother to compute 
% xsmooth = E[X(:,t) | y(:,1:T)]
% Vsmooth = Cov[X(:,t) | y(:,1:T)]
% VVsmooth = Cov[X(:,t), X(:,t-1) | y(:,1:T)] t >= 2
% loglik = sum{t=1}^T log P(y(:,t))
[xitt,xittm,Ptt,Pttm,loglik_t]=KalmanFilter(initx,initV,y',A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=KalmanSmoother(A,xitt,xittm,Ptt,Pttm,C,R);

% Create all versions of xsmooth, y, Vsmooth, VVsmooth needed
% Global + block
r = size(xsmooth);
xSmooth_C = cell(1, length(block));
for i=1:length(block)
    m = [xsmooth(1,:); xsmooth(1+i,:)];
    xSmooth_C{i} = m;
end

y_C = cell(1, length(block));
prev_split = 1;
split = 1;
for i=1:length(block)
    split = split + block(i);
    y_C{i} = y(prev_split:(split-1),:);
    prev_split = split;
end

Vsmooth_C = cell(1, length(block));
for i=1:length(block)
    V = zeros(2,2,T);
    V(1,1,:) = Vsmooth(1,1,:);
    V(2,2,:) = Vsmooth(1+i,1+i,:);
    V(1,2,:) = Vsmooth(1,1+i,:);
    V(2,1,:) = Vsmooth(1,1+i,:);
    Vsmooth_C{i} = V;
end

% All factors
xSmooth_AQ = cell(1, r(1));

for i=1:r(1)
    m = xsmooth(i,:);
    xSmooth_AQ{i} = m;
end

Vsmooth_AQ = cell(1, r(1));
for i=1:r(1)
    V = Vsmooth(i,i,:);
    Vsmooth_AQ{i} = V;
end

VVsmooth_AQ = cell(1, r(1));
for i=1:r(1)
    VV = VVsmooth(i,i,:);
    VVsmooth_AQ{i} = VV;
end

delta_C = cell(1, r(1)-1);
gamma_C = cell(1, r(1)-1);

for i=1:(length(block))
    delta = zeros(block(i), 2);
    gamma = zeros(2, 2);
    delta_C{i} = delta;
    gamma_C{i} = gamma;
end

gamma_AQ = zeros(1, r(1));
beta_AQ = zeros(1, r(1));

for t=1:T
    % C
    for i=1:length(block)
        delta_C{i} = delta_C{i} + y_C{i}(:,t)*xSmooth_C{i}(:,t)';
        gamma_C{i} = gamma_C{i} + xSmooth_C{i}(:,t)*xSmooth_C{i}(:,t)' + Vsmooth_C{i}(:,:,t);
    end
    
    % AQ
    for i=1:r(1)
        gamma_AQ(i) = gamma_AQ(i) + xSmooth_AQ{i}(:,t)*xSmooth_AQ{i}(:,t)' + Vsmooth_AQ{i}(:,:,t);
        if t>1
            beta_AQ(i) = beta_AQ(i) + xSmooth_AQ{i}(:,t)*xSmooth_AQ{i}(:,t-1)' + VVsmooth_AQ{i}(:,:,t); 
        end
    end
    
end
gamma1_AQ = zeros(1, r(1));
gamma2_AQ = zeros(1, r(1));
for i=1:r(1)
    gamma1_AQ(i) = gamma_AQ(i) - xSmooth_AQ{i}(:,T)*xSmooth_AQ{i}(:,T)' - Vsmooth_AQ{i}(:,:,T);
    gamma2_AQ(i) = gamma_AQ(i) - xSmooth_AQ{i}(:,1)*xSmooth_AQ{i}(:,1)' - Vsmooth_AQ{i}(:,:,1);
end

x1 = xsmooth(:,1);
V1 = Vsmooth(:,:,1);

