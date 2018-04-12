function [preparedData, nanMatrix, blockCount, selection] = ...
            PrepareData(data, deflate, logdiff, blockStructure)

preparedData = data;
preparedData(preparedData == 0) = NaN;

if deflate
    preparedData = deflateAdjust(preparedData);
end
% Must take YoY variables into account at some point
if logdiff
    preparedData = diff(log(preparedData));
    preparedData(preparedData == Inf | preparedData == -Inf) = NaN;
end

[~,b]=size(blockStructure);
blockCount = zeros(1,b);
selection = [];
for i=1:b
    block = blockStructure(~isnan(blockStructure(:,i)),i);
    blockCount(i) = length(block);
    selection = [selection block'];
end

preparedData = preparedData(:, selection);

prepInputMask = ~isnan(preparedData);
[T,n] = size(preparedData);
nanMatrix = zeros(n,n,T);
for t=1:T
    nanMatrix(:,:,t) = diag(prepInputMask(t,:));
end
