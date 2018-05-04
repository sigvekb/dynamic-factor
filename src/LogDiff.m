function [logdiffData] = LogDiff(logdiff, data, YoY)
logdiffData = data;
logdiffData(logdiffData == 0) = NaN;
[T, n] = size(logdiffData);

if logdiff
    diffedData = zeros(T-1,n);
    for v=1:n
        if YoY(v)
            first12 = log(logdiffData(1:(end-12),v));
            next12 = log(logdiffData(13:end,v));
            logdiff12 = next12-first12;
            diffedData(:,v) = [NaN(11,1,'like',logdiff12);logdiff12];
        else
            diffedData(:,v) = diff(log(logdiffData(:,v)));
        end
    end
    diffedData(diffedData == Inf | diffedData == -Inf) = NaN;
    logdiffData = diffedData;
end