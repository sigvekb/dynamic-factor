function [forecasts] = ForecastARIMA(data, H, ar, ma, oosm)
H_len = length(H);
maxH = max(H);

[~,n] = size(data);
forecasts = zeros(oosm,n,H_len);
ARIMAmdl = arima(ar, 0, ma);

for t=1:(oosm+maxH-1)
    fprintf('\nARIMA - forecasting month: %2d of %2d\n', t, oosm+maxH-1);
    removeMonths = oosm+maxH-t;
    forecastData_ARIMA = data(1:(end-removeMonths),:);
    
    fit = estimate(ARIMAmdl,forecastData_ARIMA,'display','off');
    [forecasted,~] = forecast(fit,maxH,'Y0',forecastData_ARIMA);

    for h=1:H_len
        horizon=H(h);
        if removeMonths>=horizon && (removeMonths-horizon)<oosm
            forecasts(t-maxH+horizon,:,h) = forecasted(horizon, :);
        end
    end
end

forecasts(forecasts == 0) = NaN;