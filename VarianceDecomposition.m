function [varDecomp] = VarianceDecomposition(normData, factors, loadings)
[n,r] = size(loadings);
varDecomp = zeros(n,4);

i_fact = normData;

all = cell(1, r);
for f=1:r
    all{f} = factors(:,f)*loadings(:,f)';
end

for v=1:n
    for f=1:r
        if f==1
            varDecomp(v,1) = var(all{f}(:,v));
            i_fact(:,v) = i_fact(:,v) - all{f}(:,v);
        else
            if sum(all{f}(:,v))~=0
                varDecomp(v,2) = var(all{f}(:,v));
                i_fact(:,v) = i_fact(:,v) - all{f}(:,v);
            end
        end
    end
    varDecomp(v,3) = nanvar(i_fact(:,v));
    varDecomp(v,4) = sum(varDecomp(v,1:3));
end