function [G, rho, Ainit] = RestrictLagMatrix(A,lags,selfLag)
Ainit = A;
r = size(lags, 2);
maxlag = max(lags);
rlag = r*maxlag;

G = zeros(rlag*r,rlag*r);

for j=1:maxlag % for each lag
    for i=1:r % each factor
        lag = lags(i);
        fRow = (j-1)*r*r+(i-1)*r;
        fCol = (j-1)*r*r+(i-1)*r;        
        if j<=lag
            if selfLag
                for v=1:r % for each factorlag
                    if i~=v
                        G(fRow+v,fCol+v) = 1;
                        Ainit(v, (j-1)*r+i) = 0;
                    end
                end
            end          
        else
            for v=1:r
            	G(fRow+v,fCol+v) = 1;
                Ainit(v, (j-1)*r+i) = 0;
            end
        end
    end
end
G( ~any(G,2), : ) = [];  % Remove all-zero rows
rho = zeros(size(G,1),1);


