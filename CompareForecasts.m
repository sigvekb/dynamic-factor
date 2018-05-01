function [relRMSE, p] = CompareForecasts(actual, forecast1, forecast2, H)
H_len = length(horizons);
relRMSE = zeros(H_len,1);
for i=1:H_len
    h = H(i);
    relRMSE(i) = RMSFE_DFM(salmonIndex,i) / RMSFE_VAR(salmonIndex,i);

    forecast_VAR = VAR_forecasts(h:end,1,i);
    error_DFM = actual_DFM(h:end,salmonIndex)-mainForecast_DFM;
    error_VAR = actual_VAR(h:end,salmonIndex)-forecast_VAR;
    [~,p] = DieboldMariano(error_DFM, error_VAR);
    DM(i) = p;
end
