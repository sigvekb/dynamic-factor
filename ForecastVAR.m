function [forecasts] = ForecastVAR(data, H, p, oosm)
H_len = length(H);
maxH = max(H);

[~,n] = size(data);
forecasts = zeros(oosm,n,H_len);
VARmdl = varm(n, p);

for t=1:(oosm+maxH-1)
    fprintf('\nVAR - forecasting month: %2d of %2d\n', t, oosm+maxH-1);
    removeMonths = oosm+maxH-t;
    forecastData_VAR = data(1:(end-removeMonths),:);
    
    fit = estimate(VARmdl,forecastData_VAR,'display','off');
    forecasted = forecast(fit,maxH,forecastData_VAR);

    for h=1:H_len
        horizon=H(h);
        if removeMonths>=horizon && (removeMonths-horizon)<oosm
            forecasts(t-maxH+horizon,:,h) = forecasted(horizon, :);
        end
    end
end

forecasts(forecasts == 0) = NaN;