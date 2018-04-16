function [x, F_hat, iter, C, A, Q] = ...
    DynamicFactorModel(X,q,r,p,max_iter, thresh, block, W)
% Extracts the unobservable factors using QML 
% Max Likelihood estimates using the Expectation Maximization (EM) algorithm 
% (Doz, Giannone and Reichlin, 2012) 
%                         
% INPUTS
% X - matrix of observable variables
% r - # of static factors
% q - # of dynamic factors
% p - # length of ar filter on common factors
% max_iter - max # of iterations in estimation
%
% OUTPUTS
% F_hat -   factors from QML

OPTS.disp = 0;
[~,n] = size(X);

demean = bsxfun(@minus, X, nanmean(X));
x = bsxfun(@rdivide, demean, nanstd(X));

% Number of dynamic factors must not exceed total number of factors
if r < q
    error('q has to be less or equal to r')
end

nlag = p-1;

A_temp = zeros(r,r*(nlag + 1))';
I = eye(r*(nlag+1),r*(nlag + 1));
A = [A_temp';I(1:end-r,1:end)];

Q = zeros((nlag+1)*r,(nlag+1)*r);
Q(1:r,1:r) = eye(r);

OPTS.disp=0;

% Create pseudo dataset, replacing NaN with 0, in order to find PCs
x_noNaN = x;
x_noNaN(isnan(x_noNaN)) = 0;

% Extract the first r eigenvectors and eigenvalues from cov(x)
[ v, ~ ] = eigs(cov(x_noNaN),r,'lm',OPTS);

chi = x_noNaN*(v*v');                       % Common component

F = x_noNaN*v;

if p > 0    
    z = F;
    Z = [];
    for kk = 1:p
        Z = [Z z(p-kk+1:end-kk,:)]; % stacked regressors (lagged SPC)
    end
    z = z(p+1:end,:);
    % run the var chi(t) = A*chi(t-1) + e(t);
    A_temp = ((Z'*Z)\Z')*z;        % OLS estimator of the VAR transition matrix
    A(1:r,1:r*p) = A_temp';
    e = z  - Z*A_temp;              % VAR residuals
    H = cov(e);                     % VAR covariance matrix
    Q(1:r,1:r) = H;                 % variance of the VAR shock when s=0
end

Q = diag(diag(Q));
A = diag(diag(A));

R = diag(diag(nancov(x-chi)));         % R diagonal

z = F;
Z = [];
for kk = 0:nlag
    Z = [Z z(nlag-kk+1:end-kk,:)];  % stacked regressors (lagged SPC)
end

% Initial factor values, initial state covariance
initx = Z(1,:)';                                                                                
initV = reshape(pinv(eye((r*nlag+1))^2-kron(A,A))*Q(:),r*(nlag+1),r*(nlag+1));

C = [v zeros(n,r*(nlag))];

% Restrict C matrix to block structure
keep = 0;
for elem=1:(length(block)-1)
    keep = keep + block(elem);
    dimV = size(v);
    for col=(elem+1):dimV(2)
        if col == elem+1
            for row = (keep+1):dimV(1)
                C(row, col) = 0;
            end
        else
            for row = 1:keep
                C(row, col) = 0;
            end
        end
    end
end

f=1; % Number of global factors
[H, K] = RestrictLoadingMatrix(n,r,f,block);

% Initialize the estimation and define ueful quantities
previous_loglik = -inf;
iter = 1;
LL = zeros(1, max_iter);
converged = 0;

y = x';

% Estimation with the EM algorithm
% Repeat the algorithm until convergence
while (iter < max_iter) && ~converged
    
    [A, C, Q, R, x1, V1, loglik, xsmooth] = ...
        EMiteration(y, A, C, Q, R, initx, initV, H, K, W);
    P1sum = V1 + x1*x1';
    % E-STEP
    % Compute the expected sufficient statistics 
    % First iteration uses initial values for A, C, Q, R, initx, initV 
    % obtained from the principal component analysis
    % 
    % beta   = sum_t=1^T (f_t * f'_t-1)
    % gamma  = sum_t=1^T (f_t * f'_t)
    % delta  = sum_t=1^T (x_t * f'_t)
    % gamma1 = sum_t=2^T (f_t-1 * f'_t-1)
    % gamma2 = sum_t=2^T (f_t * f'_t)
    % P1sum    variance of the initial state
    % x1sum    expected value of the initial state
%     [beta_AQ, delta_C, gamma1_AQ, gamma2_AQ, x1sum, V1, loglik, ~, delta, gammaKronW, gamma] = ...
%         Estep(y, A, C, Q, R, initx, initV, block, W);
                                              
    % M-STEP 
    % Compute the parameters of the model as a function of the sufficient
    % statistics
    %
    % The parameters are found through maximum likelihood optimisation
    % In the EM algorithm we substitute the sufficient statistics 
    % calculated earlier.
%     [A, C, Q, R] = ...
%         Mstep(x, p, r, block, ...
%               beta_AQ, delta_C, gamma1_AQ, gamma2_AQ, ...
%               delta, gamma, gammaKronW, R, H, K);
 
    % Update the log likelihood                                                                          
    LL(iter) = loglik;
    
    %Check for convergence
    converged = em_converged(loglik, previous_loglik, thresh,1);
    if mod(iter, 100)==0
        fprintf('Iteration %d: %f\n', iter, loglik);
    end
    
    % Set up for next iteration
    previous_loglik = loglik;
    initx = x1;
    initV = (P1sum - initx*initx');
    iter =  iter + 1;
end

[xitt,xittm,Ptt,Pttm,~]=KalmanFilter(initx,initV,x,A,C,R,Q);
[xsmooth, ~, ~]=KalmanSmoother(A,xitt,xittm,Ptt,Pttm,C,R);

F_hat =  xsmooth(1:r,:)';


function [converged, decrease] = em_converged(loglik, previous_loglik, threshold, check_increased)
% EM_CONVERGED Has EM converged?
% The EM optimisation converges if the slope of the log-likelihood 
% function falls below the set threshold, 
% i.e., |f(t) - f(t-1)| / avg < threshold,
% where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
% 'threshold' defaults to 1e-4.
%
% This stopping criterion is from Numerical Recipes in C p423

if nargin < 3, threshold = 1e-4; end
if nargin < 4, check_increased = 1; end

converged = 0;
decrease = 0;

if check_increased
    if loglik - previous_loglik < -1e-3 % allow for a little imprecision
             fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
        decrease = 1;
    end
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end
