function [H, k, Cinit] = RestrictLoadingMatrix(n, r, f, blockStruct, C)
Cinit = C;

H = zeros(n*r,n*r);
[~,blocks] = size(blockStruct);

for i=1:blocks
    block = blockStruct(:,i);
    fRow = n*i;
    fCol = n*i; 
    
    % Add restrictions
    for v=1:n    
        if ~ismember(v, block)
            H(fRow+v,fCol+v) = 1;
            Cinit(v,i+f) = 0;
        end
    end
end

H( ~any(H,2), : ) = [];  % Remove all-zero rows
k = zeros(size(H,1),1);


