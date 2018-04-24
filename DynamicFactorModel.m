function [x, F_hat, iter, Cout, Aout, Qout, R] = ...
    DynamicFactorModel(X, r, g, max_iter, thresh, blockStruct, W, VARlags, selfLag, restrictQ)
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
maxlag = max(VARlags);
rlag = r*maxlag;
VARlags = [ones(1,g-1)*VARlags(1) VARlags];

demean = bsxfun(@minus, X, nanmean(X));
x = bsxfun(@rdivide, demean, nanstd(X));

A_temp = zeros(r,rlag)';
I = eye(rlag,rlag);
A = [A_temp';I(1:end-r,1:end)];
Q = zeros(rlag,rlag);
Q(1:r,1:r) = eye(r);

% Create pseudo dataset, replacing NaN with 0, in order to find PCs
rng('default'); % Set specified seed, so that eigs is not random
x_noNaN = x;
x_noNaN(isnan(x_noNaN)) = randn();

% Extract the first r eigenvectors and eigenvalues from cov(x)
[ v, ~ ] = eigs(cov(x_noNaN),r,'lm',OPTS);
chi = x_noNaN*(v*v'); % Common component
F = x_noNaN*v;
if maxlag > 0    
    z = F;
    Z = [];
    for kk = 1:maxlag
        Z = [Z z(maxlag-kk+1:end-kk,:)]; % stacked regressors (lagged SPC)
    end
    z = z(maxlag+1:end,:);
    % run the var chi(t) = A*chi(t-1) + e(t);
    A_temp = ((Z'*Z)\Z')*z;        % OLS estimator of the VAR transition matrix
    A(1:r,1:r*maxlag) = A_temp';
    e = z  - Z*A_temp;              % VAR residuals
    H = cov(e);                     % VAR covariance matrix
    Q(1:r,1:r) = H;                 % variance of the VAR shock when s=0
    if restrictQ
        Q = diag(diag(Q));
    end
end

% Initial factor values, initial state covariance
z = F;
Z = [];
lagMinus1 = maxlag-1; %dont fully understand yet
for kk = 0:lagMinus1
    Z = [Z z(lagMinus1-kk+1:end-kk,:)];  % stacked regressors (lagged SPC)
end
initx = Z(1,:)';                                                                                
initV = cov(Z);

Atemp = A(1:r,1:r*maxlag);
[G, rho, Ainit] = RestrictLagMatrix(Atemp, VARlags, selfLag);
A(1:r,1:r*maxlag) = Ainit;

C = v;
[H, kappa, Cinit] = RestrictLoadingMatrix(n,r,g,blockStruct,C);
C = [Cinit zeros(n,r*(lagMinus1))];

R = diag(diag(nancov(x-chi)));

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
        EMiteration(y, r, A, C, Q, R, initx, initV, H, kappa, G, rho, W, restrictQ);
    
    % Update the log likelihood                                                                          
    LL(iter) = loglik;
    
    %Check for convergence
    converged = em_converged(loglik, previous_loglik, thresh);
    fprintf('Iteration %3d: %6.2f\n', iter, loglik);
    
    % Set up for next iteration
    previous_loglik = loglik;
    initx = x1;
    initV = V1;
    iter =  iter + 1;
end

F_hat =  xsmooth(1:r,:)';
Cout = C(1:n,1:r);
Aout = A(1:r,1:rlag);
Qout = Q(1:r,1:r);


function [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
% The EM optimisation converges if the slope of the log-likelihood 
% function falls below the set threshold, i.e., |f(t) - f(t-1)| / avg < threshold,
% where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is loglik at iteration t.
% 'threshold' defaults to 1e-4.

if nargin < 3, threshold = 1e-4; end

converged = 0;
decrease = 0;

if loglik - previous_loglik < -1e-3 % allow for a little imprecision
    fprintf('!! Decrease');
    decrease = 1;
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end
