function [x, F_hat, iter, C, A, Q] = ...
    DynamicFactorModel(X, q, r, p, max_iter, thresh, blockCount, blockStruct, W)
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
nlag = p-1;

% Number of dynamic factors must not exceed total number of factors
if r < q
    error('q has to be less or equal to r')
end

demean = bsxfun(@minus, X, nanmean(X));
x = bsxfun(@rdivide, demean, nanstd(X));

A_temp = zeros(r,r*(nlag + 1))';
I = eye(r*(nlag+1),r*(nlag + 1));
A = [A_temp';I(1:end-r,1:end)];
Q = zeros((nlag+1)*r,(nlag+1)*r);
Q(1:r,1:r) = eye(r);

% Create pseudo dataset, replacing NaN with 0, in order to find PCs
rng('default');
x_noNaN = x;
x_noNaN(isnan(x_noNaN)) = randn();

% Extract the first r eigenvectors and eigenvalues from cov(x)
[ v, ~ ] = eigs(cov(x_noNaN),r,'lm',OPTS);
chi = x_noNaN*(v*v'); % Common component
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

% Initial factor values, initial state covariance
z = F;
Z = [];
for kk = 0:nlag
    Z = [Z z(nlag-kk+1:end-kk,:)];  % stacked regressors (lagged SPC)
end
initx = Z(1,:)';                                                                                
initV = cov(Z);

A = diag(diag(A));
C = [v zeros(n,r*(nlag))];
Q = diag(diag(Q));
R = diag(diag(nancov(x-chi)));

f=1; % Number of global factors
[H, K, Cinit] = RestrictLoadingMatrix(n,r,f,blockStruct,C);

C = Cinit;

% Initialize the estimation and define uesful quantities
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
    
    % Update the log likelihood                                                                          
    LL(iter) = loglik;
    
    %Check for convergence
    converged = em_converged(loglik, previous_loglik, thresh,1);
    fprintf('Iteration %3d: %6.2f\n', iter, loglik);
    
    % Set up for next iteration
    previous_loglik = loglik;
    initx = x1;
    initV = V1;
    iter =  iter + 1;
end

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
             fprintf(1, '!! Decrease', previous_loglik, loglik);
        decrease = 1;
    end
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end
