function [preparedData] = prepData(data, inflateData, inflate)
% Deflate data
[T,n] = size(data);
infData = zeros(T,n);
if inflate
    for t=1:T
        for v=1:n
            infData(t,v) = data(t,v)/inflateData(t);
        end
    end
end

% Log prices, take first difference
if inflate
    preparedData = diff(log(infData));
else
    preparedData = diff(log(data));
end

