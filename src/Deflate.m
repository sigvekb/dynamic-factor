function [deflatedData] = Deflate(deflate, data) 
% Add startDate later, so different length data can be deflated
deflatedData = data;
if deflate
    inflationFile = 'InflationData.xlsx';
    [globalInflation, ~] = xlsread(inflationFile, 1);

    % Deflate data
    [T,n] = size(data);
    deflated = zeros(T,n);

    for t=1:T
        for v=1:n
            deflated(t,v) = data(t,v)/globalInflation(t);
        end
    end
end



