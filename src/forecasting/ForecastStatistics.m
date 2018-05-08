function [statistics] = ForecastStatistics(actual, f_DFM, benchmarks, H, data, oOSM)
% Calculate all necessary statistics including RMSE, RMSFE, relative RMSE,
% Jarque-Bera, Ljung-Box, Engle's test for heteroscedasticity, 
% sample covariance
H_len = size(f_DFM,2);

RMSE = zeros(1,H_len);
RMSFE = zeros(1,H_len);
Relative_RMSE = zeros(3,H_len);
DM_benchmarks = zeros(3,H_len);
LB = zeros(1,H_len);
JB = zeros(1,H_len);
Engle = zeros(1,H_len);
sampleCorr  = zeros(1,H_len);
stdErr = zeros(1,H_len);
relCorr = zeros(3,H_len);
MASE = zeros(1,H_len);

error_DFM = bsxfun(@minus, f_DFM, actual);

naive_MAE = (1/(size(data,1)-oOSM-1))*...
            nansum(abs(data(2:(end-oOSM), 1)-data(1:(end-oOSM-1), 1)));

for h=1:H_len
    horizon = H(h);
    [rmsfe, rmse] = CalcRMSFE(actual, f_DFM(:,h), horizon);
    RMSE(1,h) = rmse;
    RMSFE(1,h) = rmsfe;
    s = corrcoef(f_DFM(:,h),actual);
    sampleCorr(1,h) = s(2,1);
    for i=1:3
        [~, rmse] = CalcRMSFE(actual, benchmarks(:,h,i), horizon);
        Relative_RMSE(i,h) = RMSE(1,h) / rmse; 
        errorBenchmark = bsxfun(@minus, benchmarks(:,h,i), actual);
        [~,p] = DieboldMariano(error_DFM(:,h), errorBenchmark);
        DM_benchmarks(i,h) = p;
        
        b = corrcoef(benchmarks(:,h,i),actual);
        benchCorr = b(2,1);
        relCorr(i,h) = sampleCorr(1,h) / benchCorr;
    end
    
    %Ljung-Box Q test
    [~,p_lb] = lbqtest(error_DFM(:,h));
    %Jarque Bera test
    [~,p_jb]= jbtest(error_DFM(:,h));
    %Engle test for heteroscedastisity
    [~,p_en] = archtest(error_DFM(:,h));
    %Sample correlation
    s = corrcoef(f_DFM(:,h),actual);
    sampleCorr(1,h) = s(2,1);
    % Standard Error
    stdErr(1,h) = std(error_DFM(:,h));
        
    LB(h) = p_lb;
    JB(h) = p_jb;
    Engle(h) = p_en;
    
    forecastError = (1/oOSM)*nansum(abs(error_DFM(:,h)));
    MASE(1,h) = forecastError/naive_MAE;
end

% Calculate MASE


statistics = [RMSE;
              RMSFE;
              Relative_RMSE;
              DM_benchmarks;
              LB;
              JB;
              Engle;
              sampleCorr;
              stdErr;
              relCorr;
              MASE];
