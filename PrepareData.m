function [preparedData, nanMatrix, newBlockStruct, blockCount, varsFound] = ...
            PrepareData(data, deflate, logdiff, blockStruct, YoY)

preparedData = data;
preparedData(preparedData == 0) = NaN;
[T, n] = size(preparedData);

if deflate
    preparedData = deflateAdjust(preparedData);
end
% Must take YoY variables into account at some point

if logdiff
    diffedData = zeros(T-1,n);
    for v=1:n
        if YoY(v)
            first12 = log(preparedData(1:(end-12),v));
            next12 = log(preparedData(13:end,v));
            logdiff12 = next12-first12;
            diffedData(:,v) = [NaN(11,1,'like',logdiff12);logdiff12];
        else
            diffedData(:,v) = diff(log(preparedData(:,v)));
        end
    end
    diffedData(diffedData == Inf | diffedData == -Inf) = NaN;
    preparedData = diffedData;
end

% Define block structure in terms of new ordering
[d,b]=size(blockStruct);
newBlockStruct = zeros(d,b);
varsFound = [];

blockCount = zeros(1,b);
for i=1:b
    block = blockStruct(~isnan(blockStruct(:,i)),i);
    bLen = length(block);
    for j=1:bLen
        var = blockStruct(j,i);
        found = find(varsFound==var);
        if isempty(found)
            newBlockStruct(j,i) = length(varsFound)+1;
            varsFound = [varsFound var];
        else
            newBlockStruct(j,i) = found;
        end
    end
    blockCount(i) = bLen;
end

preparedData = preparedData(:, varsFound);

prepInputMask = ~isnan(preparedData);
[T,n] = size(preparedData);
nanMatrix = zeros(n,n,T);
for t=1:T
    nanMatrix(:,:,t) = diag(prepInputMask(t,:));
end
