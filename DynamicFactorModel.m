function [x,F_hat,iter, C, A, Q] = DynamicFactorModel(X,q,r,p,max_iter, thresh, block)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "A Quasi?Maximum Likelihood Approach for Large, Approximate Dynamic Factor Models," 
% The Review of Economics and Statistics, MIT Press, vol. 94(4), pages 1014-1024, November 2012.
% Catherine Doz, Universite' Cergy-Pontoise
% Domenico Giannone, Universite' Libre de Bruxelles, ECARES and CEPR
% Lucrezia Reichlin, London Business School and CEPR 
%
%
% Programs are also available at: http://homepages.ulb.ac.be/~dgiannon/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% DynFA:             extracts the unobservable factors using three different methods 
%                         
%       - QML:      Max Likelihood estimates using the Expectation Maximization (EM) algorithm 
%                    (Doz, Giannone and Reichlin, 2012) 
%                         
%       - TWO STEP: Principal components + Kalman filtering 
%                   Doz, Catherine & Giannone, Domenico & Reichlin, Lucrezia, 2011.
%                   "A two-step estimator for large approximate dynamic factor models based on Kalman filtering," 
%                   Journal of Econometrics, Elsevier, vol. 164(1), pages 188-205, September.
% 
%       - PC:       principal components 
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
% F_pc  -   factors using principal components
% F_kal -   factors using from two steps

OPTS.disp = 0;
[T,N] = size(X);

demean = bsxfun(@minus, X, nanmean(X));
x = bsxfun(@rdivide, demean, nanstd(X));


% the number of static factors cannot be greater than dynamic
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

% Create pseudo dataset, replacing NaN with 0, in order to make
% the initialization good
pseudo_x = x;
pseudo_x(isnan(pseudo_x)) = 0;

% Extract the first r eigenvectors and eigenvalues from cov(x)
[ v, ~ ] = eigs(cov(pseudo_x),r,'lm',OPTS);

chi = pseudo_x*(v*v');                       % Common component

F = pseudo_x*v;

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

C = [v zeros(N,r*(nlag))];

% Restric C matrix to block structure
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

% Initialize the estimation and define ueful quantities
previous_loglik = -inf;
iter = 1;
LL = zeros(1, max_iter);
converged = 0;

y = x';

% Estimation with the EM algorithm
% Repeat the algorithm until convergence
while (iter < max_iter) && ~converged
    
    % E step : computes the expected sufficient statistics 
    % First iteration uses initial values for A, C, Q, R, initx, initV 
    % obtained from the principal component analysis
    [beta_t_AQ, gamma_t_C, delta_t_C, gamma1_t_AQ, gamma2_t_AQ, x1, V1, loglik_t,xsmooth] = ...
        Estep(y, A, C, Q, R, initx, initV, block);
    
    % Hand over the expected sufficient statistics from E-step
    beta_AQ = beta_t_AQ;        % beta   = sum_t=1^T (f_t * f'_t-1)
    gamma_C = gamma_t_C;        % gamma  = sum_t=1^T (f_t * f'_t)   
    delta_C = delta_t_C;        % delta  = sum_t=1^T (x_t * f'_t)
    gamma1_AQ = gamma1_t_AQ;    % gamma1 = sum_t=2^T (f_t-1 * f'_t-1)
    gamma2_AQ = gamma2_t_AQ;    % gamma2 = sum_t=2^T (f_t * f'_t)
    P1sum = V1 + x1*x1';        % P1sum    variance of the initial state
    x1sum = x1;                 % x1sum    expected value of the initial state
                                                                            
    % Update the log likelihood                                                                           
    loglik = loglik_t;
    LL(iter) = loglik;
    
    %%% M step 
    % Compute the parameters of the model as a function of the sufficient
    % statistics
    
    % The parameters are found through maximum likelihood optimisation
    % In the EM algorithm we substitute the sufficient statistics 
    % calculated earlier.
    
    % C = (sum_t=1^T x_t*f'_t)* (sum_t=1^T f_t*f'_t)^-1 
    prev_split = 1;
    split = 1;
    for i=1:length(block)
        split = split + block(i);
        C_temp = delta_C{i} * pinv(gamma_C{i});
        C(prev_split:(split-1), 1) = C_temp(:,1);
        C(prev_split:(split-1), 1+i) = C_temp(:,2);
        prev_split = split;
    end
    
    if p > 0
    
        % A = (sum_t=2^T f_t*f'_t-1)* (sum_2=1^T f_t-1*f'_t-1)^-1
        Atemp_AQ = zeros(1, r);
        for i=1:r
            Atemp_AQ(i) = beta_AQ(i) / gamma1_AQ(i);
            A(i,i) = Atemp_AQ(i);
        end
        
        % Q = ( (sum_t=2^T f_t*f'_t) - A * (sum_2=1^T f_t-1*f'_t) )/(T-1)
        for i=1:r
            H = (gamma2_AQ(i)-Atemp_AQ(i)*beta_AQ(i)) / (T-1);
            Q(i,i) = H;
        end
    end
    
    % R = ( sum_t=1^T (x_t*x'_t) - C * f_t*x'_t) )/T 
    prev_split = 1;
    split = 1;
    Cdim = size(C);
    delta = zeros(Cdim(1), Cdim(2));
    for i=1:length(block)
        split = split + block(i);
        delta(prev_split:(split-1), 1) = delta_C{i}(:,1);
        delta(prev_split:(split-1), 1+i) = delta_C{i}(:,2);
        prev_split = split;
    end
    R = (x'*x - C*delta')/T;
    
    RR = diag(R); RR(RR<1e-7) = 1e-7; R = diag(RR); 
      
    R = diag(diag(R));                  % R diagonal
    
    if mod(iter, 100)==0
        fprintf('Iteration %d: %f\n', iter, loglik);
    end
    
    initx = x1sum;
    initV = (P1sum - initx*initx');
    
    converged = em_converged(loglik, previous_loglik, thresh,1);
    
    iter =  iter + 1;
    
    previous_loglik = loglik;
end

[xitt,xittm,Ptt,Pttm,loglik_t]=KalmanFilter(initx,initV,x,A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=KalmanSmoother(A,xitt,xittm,Ptt,Pttm,C,R);

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
