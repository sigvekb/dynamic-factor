function [varDecomp] = VarianceDecomposition(normData, factors, loadings)
[n,r] = size(loadings);
varDecomp = zeros(n,r+2);

% Find variance explained for each factor
for f=1:r
    factorVar = var(factors(:,f));
    varDecomp(:,f) = loadings(:,f).^2 * factorVar;
end

% Sum over variances of block factors
varDecomp(:,r+1) = sum(varDecomp(:,2:r),2);

% Find error term variance
idio = normData - factors*loadings';
varDecomp(:,r+2) = nanvar(idio)';