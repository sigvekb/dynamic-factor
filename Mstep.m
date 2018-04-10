function [A, C, Q, R] = Mstep(x, p, r, block, beta_AQ, gamma_C, delta_C, ...
                              gamma1_AQ, gamma2_AQ)
                          

[T, ~] = size(x);
% Update C (loadings matrix)
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

% If we use dynamic factors, we update the prediction equation parameters
if p > 0
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