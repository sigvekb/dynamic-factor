function [varDecomp] = VarianceDecomposition(normData, factors, loadings, g)
[n,r] = size(loadings);
varDecomp = zeros(n,r+3);

% Find variance explained for each factor
for f=1:r
    factorVar = var(factors(:,f));
    varDecomp(:,f) = loadings(:,f).^2 * factorVar;
end

% Create summary columns
% Sum over variances of global factors
varDecomp(:,r+1) = sum(varDecomp(:,1:g),2);

% Sum over variances of block factors
varDecomp(:,r+2) = sum(varDecomp(:,(g+1):r),2);

% Find error term variance
idio = normData - factors*loadings';
varDecomp(:,r+3) = nanvar(idio)';