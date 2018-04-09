function [deflated] = deflateAdjust(data)
inflationFile = "InflationData.xlsx";
[globalInflation, ~] = xlsread(inflationFile, 1);

% Deflate data
[T,n] = size(data);
deflated = zeros(T,n);

for t=1:T
    for v=1:n
        deflated(t,v) = data(t,v)/globalInflation(t);
    end
end

