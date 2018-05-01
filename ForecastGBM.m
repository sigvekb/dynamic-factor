function [forecasts] = ForecastGBM(data, H, oosm)
H_len = length(H);
maxH = max(H);

forecasts = zeros(oosm,H_len);


fprintf('\nGBM - forecasting');
for t=1:(oosm+maxH-1)
    removeMonths = oosm+maxH-t;
    forecastData = data(1:(end-removeMonths),:);
    
    mu = nanmean(forecastData);
    sigma = nanstd(forecastData);

    for h=1:H_len
        horizon=H(h);
        if removeMonths>=horizon && (removeMonths-horizon)<oosm
            forecasts(t-maxH+horizon,h) = randn(1,1)*sigma + mu*horizon;
        end
    end
end