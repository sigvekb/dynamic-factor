function [varDecomp] = DynFactorTest(normData, factors, loadings)
[T,n] = size(normData);
[~,r] = size(loadings);
varDecomp = zeros(n,4);

i_fact = zeros(T,n);
i_fact = i_fact + normData;

all = cell(1, r);
all{1} = factors(:,1)*loadings(:,1)';
for f=2:r
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
    varDecomp(v,3) = var(i_fact(:,v));
    varDecomp(v,4) = sum(varDecomp(v,1:3));
end