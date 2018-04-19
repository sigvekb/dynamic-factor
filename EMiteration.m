function [A, C, Q, R, x1, V1, loglik_t, xsmooth] = ...
    EMiteration(y, A, C, Q, R, initx, initV, H, K, W)

[n, T] = size(y);

[xitt,xittm,Ptt,Pttm,loglik_t, Kgain]=KalmanFilter(initx,initV,y',A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=KalmanSmoother(A,xitt,xittm,Kgain,Ptt,Pttm,C);

x1 = xsmooth(:,1);
V1 = Vsmooth(:,:,1);

% Create all versions of xsmooth, y, Vsmooth, VVsmooth needed
[r, ~] = size(xsmooth);

xSmooth_AQ = cell(1, r);
Vsmooth_AQ = cell(1, r);
VVsmooth_AQ = cell(1, r);
gamma_AQ = zeros(1, r);
beta_AQ = zeros(1, r);
for i=1:r
    xSmooth_AQ{i} = xsmooth(i,:);
    Vsmooth_AQ{i} = Vsmooth(i,i,:);
    VVsmooth_AQ{i} = VVsmooth(i,i,:);
end

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
gammaKronW = zeros(n*r, n*r);
beta = zeros(r, r);
y_noNaN = y;
y_noNaN(isnan(y_noNaN)) = 0;
for t=1:T
    deltaChange = W(:,:,t)*y_noNaN(:,t)*xsmooth(:,t)';
    delta = delta + deltaChange; % D
    
    gammaChange = xsmooth(:,t)*xsmooth(:,t)' + Vsmooth(:,:,t);
    gammaKronW = gammaKronW + kron(gammaChange, W(:,:,t));
    
    gamma = gamma + gammaChange;
    
    if t>1
        beta = beta + xsmooth(:,t)*xsmooth(:,t-1)' + VVsmooth(:,:,t); % 
    end
end

gamma1 = gamma - xsmooth(:,T)*xsmooth(:,T)' - Vsmooth(:,:,T); % 
gamma2 = gamma - xsmooth(:,1)*xsmooth(:,1)' - Vsmooth(:,:,1);

% Start M-step
%-------------------------------------------------------------------------                        
% Update C (See Banbura(2010) and Bork(2008))
% C is block-restricted and handles missing values
gammaKronR = kron(inv(gamma), R);
Cvec = gammaKronW \ delta(:);
Cres = Cvec + ((gammaKronR * H') / (H * gammaKronR * H')) * (K - H * Cvec);
C = reshape(Cres, [n,r]);
C(abs(C)<1e-10) = 0; % Remove almost-zero entries..

% Update R (See Banbura(2010))
[first, andre, tredje, fjerde, femte] = deal(zeros(n,n));
I_n = eye(n);
for t=1:T
    select = W(:,:,t);
    x_t = xsmooth(:,t);
    y_t = y_noNaN(:,t);
    
    first = first + select * (y_t * y_t') * select';
    andre = andre + select * y_t * x_t' * C' * select;
    tredje = tredje + select * C * x_t * y_t' * select;
    gamma_t = x_t * x_t' + Vsmooth(:,:,t);
    fjerde = fjerde + select*C*gamma_t*C'*select;
    femte = femte + (I_n - select) * R * (I_n - select);
end
R = 1/T * (first - andre - tredje + fjerde + femte);
R = diag(diag(R));

% NB! The A matrix is only square if the factors are all AR(1) processes
% This must be updated otherwise!
A = zeros(r,r);
Q = zeros(r,r);

% Update A (Factor VAR coefficients)
% A = (sum_t=2^T f_t*f'_t-1)* (sum_2=1^T f_t-1*f'_t-1)^-1
Atemp_AQ = zeros(1, r);
for i=1:r
    Atemp_AQ(i) = beta_AQ(i) / gamma1_AQ(i);
    A(i,i) = Atemp_AQ(i);
end

% Update Q (Factor VAR error covariance matrix)
% Q = ( (sum_t=2^T f_t*f'_t) - A * (sum_2=1^T f_t-1*f'_t) )/(T-1)
for i=1:r
    H = (gamma2_AQ(i)-Atemp_AQ(i)*beta_AQ(i)) / (T-1);
    Q(i,i) = H;
end