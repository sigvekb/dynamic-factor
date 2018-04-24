function [nanMatrix] = CreateNaNMatrix(data)
prepInputMask = ~isnan(data);
[T,n] = size(data);
nanMatrix = zeros(n,n,T);
for t=1:T
    nanMatrix(:,:,t) = diag(prepInputMask(t,:));
end