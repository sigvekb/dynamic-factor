function [MASE] = CalcMASE(r_in,rhat_error)
% Calculate MASE for returns, with naive no-change (ARIMA(0,1,0)) for price    
naive_MAE = (1/size(r_in,1))*nansum(abs(r_in));

forecastError = (1/size(rhat_error,1))*nansum(abs(rhat_error));
MASE = forecastError/naive_MAE;