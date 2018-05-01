function [forecasts, varDecomp] = ForecastDFM(data,H,oosm,...
                               g,iter,thresh,sLag,resQ,blockStruct,lags)
H_len = length(H);
maxH = max(H);

[~,n] = size(data);
forecasts = zeros(oosm,n,H_len);

for t=1:(oosm+maxH-1)
    fprintf('\nDFM - forecasting month: %2d of %2d\n', t, oosm+maxH-1);
    removeMonths = oosm+maxH-t;
    
    forecastData = [data(1:(end-removeMonths),:); NaN([maxH,n])];
    [~,factors,~,~,C,~,~,~] = ...
            DynamicFactorModel(forecastData, g, iter, thresh, sLag, resQ, blockStruct, lags, []);

    for h=1:H_len
        horizon=H(h);
        if removeMonths>=horizon && (removeMonths-horizon)<oosm
            f_index = size(factors,1)-maxH+horizon;
            forecasts(t-maxH+horizon,:,h) = (C*factors(f_index,:)')';
        end
    end 
end

% Variance Decomposition of DFM
demean = bsxfun(@minus, data, nanmean(data));
DFMnorm = bsxfun(@rdivide, demean, nanstd(demean));  
[varDecomp] = VarianceDecomposition(DFMnorm, factors(1:(end-maxH+1),:), C, g);

forecasts(forecasts == 0) = NaN;