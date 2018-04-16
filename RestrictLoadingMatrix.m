function [H, k] = RestrictLoadingMatrix(n, r, f, block)

k = zeros(n*(r-(f+1)),1);

H = zeros(n*(r-(f+1)), n*r);
b = length(block);
sumB = 0;
for i=1:b
    len = block(i);
    
    
    cols = (1:n) + (f+i-1)*n;
    zerocols = (1:len)+sumB;
    cols(zerocols) = [];
    rows = (1:(n-len)) + n*(i-1)-sumB;
    sumB = sum(block(1:i));
    
    for j = 1:length(rows)
        H(rows(j), cols(j)) = 1;
    end
end


