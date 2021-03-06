function [A, C, Q, R, x1, V1, loglik_t, xsmooth] = ...
    EMiteration(y, r, A, C, Q, R, initx, initV, H, kappa, G, rho, W, restrictQ)

[n, T] = size(y);
rlag = size(A,1);

[xitt,xittm,Ptt,Pttm,loglik_t, Kgain]=KalmanFilter(initx,initV,y',A,C,R,Q);
[xsmooth, Vsmooth, VVsmooth]=KalmanSmoother(A,xitt,xittm,Kgain,Ptt,Pttm,C);

x1 = xsmooth(:,1);
V1 = Vsmooth(:,:,1);

% Calculate necessary moments
delta = zeros(n,rlag);
gamma = zeros(rlag,rlag);
gammaKronW = zeros(n*r, n*r);
beta = zeros(rlag, rlag);
y_noNaN = y;
y_noNaN(isnan(y_noNaN)) = 0;
for t=1:T
    deltaChange = W(:,:,t)*y_noNaN(:,t)*xsmooth(:,t)';
    delta = delta + deltaChange;
    
    gammaChange = xsmooth(:,t)*xsmooth(:,t)' + Vsmooth(:,:,t);
    gamma = gamma + gammaChange;
    gammaKronW = gammaKronW + kron(gammaChange(1:r,1:r), W(:,:,t));
    
    if t>1
        beta = beta + xsmooth(:,t)*xsmooth(:,t-1)' + VVsmooth(:,:,t);
    end
end
gamma1 = gamma - xsmooth(:,T)*xsmooth(:,T)' - Vsmooth(:,:,T);
gamma2 = gamma - xsmooth(:,1)*xsmooth(:,1)' - Vsmooth(:,:,1);

% Start M-step
%-------------------------------------------------------------------------                        
% Update C (See Banbura(2010) and Bork(2008))
% C is block-restricted and handles missing values
gammaInv = CholeskyInversion(gamma(1:r,1:r));
gammaKronR = kron(gammaInv, R);

deltaC = delta(:,1:r);
gKW_Inv = CholeskyInversion(gammaKronW);
Cvec = gKW_Inv * deltaC(:);

HgH_inv = CholeskyInversion(H*gammaKronR*H');
Cres = Cvec + gammaKronR * H' * HgH_inv * (kappa-H*Cvec);

C(1:n,1:r) = reshape(Cres, [n,r]);
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

% Update A (Factor VAR coefficients)
% A = (sum_t=2^T f_t*f'_t-1)* (sum_2=1^T f_t-1*f'_t-1)^-1
gamma1Inv = CholeskyInversion(gamma1(1:rlag,1:rlag));
Atemp = beta(1:r,1:rlag) * gamma1Inv;

gamma1KronQ = kron(gamma1Inv, Q(1:r,1:r));

GgG_inv = CholeskyInversion(G * gamma1KronQ * G');
%GgG = G * gamma1KronQ * G';
Avec = Atemp(:) + gamma1KronQ * G' * GgG_inv * (rho-G*Atemp(:));
Avec(abs(Avec)<1e-10) = 0; % Remove almost-zero entries..
A(1:r,1:rlag) = reshape(Avec, [r,rlag]);

% Update Q (Factor VAR error covariance matrix)
% Q = ( (sum_t=2^T f_t*f'_t) - A * (sum_2=1^T f_t-1*f'_t) )/(T-1)
Qupdate = (gamma2(1:r,1:r) - Atemp*beta(1:r,1:rlag)') / (T-1);
if restrictQ
    Q(1:r,1:r) = diag(diag(Qupdate));
else
    Q(1:r,1:r) = Qupdate;
end