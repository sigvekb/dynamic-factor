function [x,F_hat,F_pc,F_kal,num_iter, C, A, Q] = DynFactorDoz2(X,q,r,p,max_iter, block)
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


thresh = 1e-7;
OPTS.disp = 0;
[T,N] = size(X);

Mx = mean(X);
Wx = (std(X));
% x is the matrix of standardized observable variables
x = (X-kron(ones(T,1),Mx))*diag(1./Wx);

% the number of static factors cannot be great of the dynamic ones
if r < q
    error('q has to be less or equal to r')
end

nlag = p-1;

A_temp = zeros(r,r*(nlag + 1))';
I = eye(r*(nlag+1),r*(nlag+1));
A = [A_temp';I(1:end-r,1:end)];

Q = zeros((nlag+1)*r,(nlag+1)*r);
Q(1:r,1:r) = eye(r);

OPTS.disp=0;

%extract the first r eigenvectors and eigenvalues from cov(x)
[ v, ~ ] = eigs(cov(x),r,'lm',OPTS);

chi = x*(v*v');                       % common component

d = eye(r);

F = x*v;

F_pc = F;                           % factors using principal components

if p > 0    
    z = F;
    Z = [];
    for kk = 1:p
        Z = [Z z(p-kk+1:end-kk,:)]; % stacked regressors (lagged SPC)
    end
    z = z(p+1:end,:);
    % run the var chi(t) = A*chi(t-1) + e(t);
    A_temp = inv(Z'*Z)*Z'*z;        % OLS estimator of the VAR transition matrix
    A(1:r,1:r*p) = A_temp';
    e = z  - Z*A_temp;              % VAR residuals
    H = cov(e);                     % VAR covariance matrix
    Q(1:r,1:r) = H;                 % variance of the VAR shock when s=0
end

Q = diag(diag(Q));
A = diag(diag(A));

R = diag(diag(cov(x-chi)));         % R diagonal

z = F;
Z = [];
for kk = 0:nlag
    Z = [Z z(nlag-kk+1:end-kk,:)];  % stacked regressors (lagged SPC)
end

initx = Z(1,:)';                    %initial state mean                                                            


% initial state covariance
%initV = reshape(pinv(eye((r*nlag+1))^2-kron(A,A))*Q(:),r*(nlag+1),r*(nlag+1));eye(r*(nlag+1));
initV = cov(Z);


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

% initialize the estimation and define ueful quantities

previous_loglik = -inf;

loglik = 0;
num_iter = 0;
LL = -inf;

os = size(C,1);     % number of cross sections ( N )
ss = size(A,1);     % number of factors ( r )
y = x';



converged = 0;

% estimation of the factors with the Kalman filter using as initial values
% for A, C, Q, R, initx, initV the ones computed with the principal
% components
[xitt,xittm,Ptt,Pttm,loglik_t]=K_filter(initx,initV,x,A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=K_smoother(A,xitt,xittm,Ptt,Pttm,C,R);

F_kal =  xsmooth(1:r,:)';


% factors estimation with the EM algorithm

%repeat the algorithm until convergence
while (num_iter < max_iter) && ~converged
    
    %%% E step : compute the sufficient statistics 
    
    % In our model the sufficient statistics are
    % delta = sum_t=1^T (x_t * f'_t)
    % gamma = sum_t=1^T (f_t * f'_t)   
    % gamma1 = sum_t=2^T (f_t-1 * f'_t-1)
    % gamma2 = sum_t=2^T (f_t * f'_t)
    % beta = sum_t=1^T (f_t * f'_t-1)
    % P1sum  variance of the initial state
    % x1sum  expected value of the initial state 
    
    % use the function Estep to update the expected sufficient statistics
    % note that at the first iteration  we use as initial values for A, C, Q, 
    % R, initx, initV the ones computed with the principal components
    [beta_t_AQ, gamma_t_C, gamma_t_AQ, delta_t_C, gamma1_t_AQ, gamma2_t_AQ, x1, V1, loglik_t,xsmooth] = ...
        Estep(y, A, C, Q, R, initx, initV, block);
    
    % fix the expected sufficient statistics equal to the one computed with
    % the function Estep
    beta_AQ = beta_t_AQ;        % beta   = sum_t=1^T (f_t * f'_t-1)
    gamma_C = gamma_t_C;        % gamma  = sum_t=1^T (f_t * f'_t)   
    gamma_AQ = gamma_t_AQ;      % gamma  = sum_t=1^T (f_t * f'_t)
    delta_C = delta_t_C;        % delta  = sum_t=1^T (x_t * f'_t)
    gamma1_AQ = gamma1_t_AQ;    % gamma1 = sum_t=2^T (f_t-1 * f'_t-1)
    gamma2_AQ = gamma2_t_AQ;    % gamma2 = sum_t=2^T (f_t * f'_t)
    P1sum = V1 + x1*x1';        % P1sum    variance of the initial state
    x1sum = x1;                 % x1sum    expected value of the initial state
    
                                                                            
    % update the loglikelihood                                                                           
    loglik = loglik_t;
    
    % update the counter for the iterations
    num_iter =  num_iter + 1;
    
    
    %%% M step 
    % compute the parameters of the model as a function of the sufficient
    % statistics (computed with the function Estep)
    
    % The formulas for the parameters derive from the maximum likelihood
    % method. In the EM algorithm we substitute in the ML estimator the 
    % sufficient statistics (computed in the E step) and then iterate the 
    % procedure till the maximization of the likelihood
    
    % C = (sum_t=1^T x_t*f'_t)* (sum_t=1^T f_t*f'_t)^-1 
    % substituting for the sufficient statistics
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
    loglik = loglik_t;
    
    LL = [LL loglik];
    
    fprintf('Iteration %d: %f\n', num_iter, loglik);
    
    initx = x1sum;
    initV = (P1sum - initx*initx');
    
    converged = em_converged(loglik, previous_loglik, thresh,1);
    
    previous_loglik = loglik;
    

end



[xitt,xittm,Ptt,Pttm,loglik_t]=K_filter(initx,initV,x,A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=K_smoother(A,xitt,xittm,Ptt,Pttm,C,R);

chi = xsmooth'*C'*diag(Wx) + kron(ones(T,1),Mx);

F_hat =  xsmooth(1:r,:)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [beta_AQ, gamma_C, gamma_AQ, delta_C, gamma1_AQ, gamma2_AQ, x1, V1, loglik_t, xsmooth] = ...
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

[os, T] = size(y); % os - number of input variables
ss = length(A); % Number of factors

% use the Kalman smoother to compute 
% xsmooth = E[X(:,t) | y(:,1:T)]
% Vsmooth = Cov[X(:,t) | y(:,1:T)]
% VVsmooth = Cov[X(:,t), X(:,t-1) | y(:,1:T)] t >= 2
% loglik = sum{t=1}^T log P(y(:,t))
[xitt,xittm,Ptt,Pttm,loglik_t]=K_filter(initx,initV,y',A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=K_smoother(A,xitt,xittm,Ptt,Pttm,C,R);

% Create all versions of xsmooth, y, Vsmoot, VVsmooth needed
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

% compute the expected sufficient statistics
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [converged, decrease] = em_converged(loglik, previous_loglik, threshold, check_increased)
% EM_CONVERGED Has EM converged?
% [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
%
% We have converged if the slope of the log-likelihood function falls below 'threshold', 
% i.e., |f(t) - f(t-1)| / avg < threshold,
% where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
% 'threshold' defaults to 1e-4.
%
% This stopping criterion is from Numerical Recipes in C p423
%
% If we are doing MAP estimation (using priors), the likelihood can decrase,
% even though the mode of the posterior is increasing.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xitt,xittm,Ptt,Pttm,loglik]=K_filter(initx,initV,x,A,C,R,Q)
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
% Pttm = Cov[X(:,t) | y(:,1:t-1)]
% xitt = E[X(:,t) | y(:,1:t)]
% Ptt = Cov[X(:,t) | y(:,1:t)]
%loglik - value of the loglikelihood

[T,N]=size(x);
r=size(A,1);

y=x';


xittm=[initx zeros(r,T)];
xitt=zeros(r,T);

Pttm=zeros(r,r,T);
Pttm(:,:,1)=initV;
Ptt=zeros(r,r,T);

for j=1:T
    
    L=inv(C*Pttm(:,:,j)*C'+R);
    nv=size(Pttm,1);
    xitt(:,j)=xittm(:,j)+Pttm(:,:,j)*C'*L*(y(:,j)-C*xittm(:,j));
    Ptt(:,:,j)=Pttm(:,:,j)-Pttm(:,:,j)*C'*L*C*Pttm(:,:,j); 
    xittm(:,j+1)=A*xitt(:,j);
    Pttm(:,:,j+1)=A*Ptt(:,:,j)*A'+Q;
    lik(j)=((2*pi)^(-N/2))*(abs((det(C*Pttm(:,:,j)*C'+R)))^(-.5))*...
        exp(-1/2*(y(:,j)-C*xittm(:,j))'*L*(-1/2*(y(:,j)-C*xittm(:,j))));
    
    
    e = y(:,j) - C*xittm(:,j); % error (innovation)
    n = length(e);
    ss = length(A);
    d = size(e,1);
    S = C*Pttm(:,:,j)*C' + R;
    GG = C'*diag(1./diag(R))*C;
    Sinv = inv(S);
    
    %%%%%%%%%%%%%%%%%%%%%%
    
    detS = prod(diag(R))*det(eye(ss)+Pttm(:,:,j)*GG);
    denom = (2*pi)^(d/2)*sqrt(abs(detS));
    mahal = sum(e'*Sinv*e,2);
    logl(j) = -0.5*mahal - log(denom);
    
end

loglik=sum(logl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xitT,PtT,PtTm]=K_smoother(A,xitt,xittm,Ptt,Pttm,C,R)
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
    J(:,:,i)=Ptt(:,:,i)*A'*inv(Pttm(:,:,i+1));
end

for i=1:T
    L(:,:,i)=inv(C*Pttm(:,:,i)*C'+R);
    K(:,:,i)=Pttm(:,:,i)*C'*L(:,:,i);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XC = center(X)
%CENTER XC = center(X)
%	Centers each column of X.

%	J. Rodrigues 26/IV/97, jrodrig@ulb.ac.be

[T, n] = size(X);
XC = X - ones(T,1)*(sum(X)/T); % Much faster than MEAN with a FOR loop
