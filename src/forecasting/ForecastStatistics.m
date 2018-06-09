function [statistics] = ForecastStatistics(actual, f_DFM, benchmarks, H, data, oOSM)
% Calculate all necessary statistics including RMSE, RMSFE, relative RMSE,
% Jarque-Bera, Ljung-Box, Engle's test for heteroscedasticity, 
% sample covariance
H_len = size(f_DFM,2);

RMSE = zeros(4,H_len);
Relative_RMSE = zeros(3,H_len);
DM_benchmarks = zeros(3,H_len);
LB = zeros(4,H_len);
JB = zeros(4,H_len);
Engle = zeros(4,H_len);
sampleCorr  = zeros(4,H_len);
stdErr = zeros(1,H_len);
MASE = zeros(4,H_len);
HitRate = zeros(4,H_len);

error_DFM = bsxfun(@minus, f_DFM, actual);

for h=1:H_len
    horizon = H(h);
    [~, rmse] = CalcRMSFE(actual, f_DFM(:,h), horizon);
    RMSE(1,h) = rmse;
    
    [~,p_lb] = lbqtest(error_DFM(:,h)); 
    LB(1,h) = p_lb;
    [~,p_jb]= jbtest(error_DFM(:,h));
    JB(1,h) = p_jb;
    [~,p_en] = archtest(error_DFM(:,h));
    Engle(1,h) = p_en;
    stdErr(1,h) = std(error_DFM(:,h));
    
    s = corrcoef(f_DFM(:,h),actual);
    sampleCorr(1,h) = s(2,1);
    
    MASE(1,h) = CalcMASE(data(1:(end-oOSM),1), error_DFM(:,h));
    HitRate(1,h) = CalcHR(actual,f_DFM(:,h));
    
    for i=1:3
        [~, rmse] = CalcRMSFE(actual, benchmarks(:,h,i), horizon);
        RMSE(i+1,h) = rmse;
        Relative_RMSE(i,h) = RMSE(1,h) / rmse; 
        errorBenchmark = benchmarks(:,h,i) - actual;
        [~,p] = DieboldMariano(error_DFM(:,h), errorBenchmark);
        DM_benchmarks(i,h) = p;
        
        [~,p_lb] = lbqtest(errorBenchmark); 
        LB(i+1,h) = p_lb;
        [~,p_jb]= jbtest(errorBenchmark);
        JB(i+1,h) = p_jb;
        [~,p_en] = archtest(errorBenchmark);
        Engle(i+1,h) = p_en;
        
        b = corrcoef(benchmarks(:,h,i),actual);
        sampleCorr(i+1,h) = b(2,1);
        
        HitRate(i+1,h) = CalcHR(actual,benchmarks(:,h,i));
        MASE(i+1,h) = CalcMASE(data(1:(end-oOSM),1), errorBenchmark);
    end
end
statistics = [RMSE;
              Relative_RMSE;
              DM_benchmarks;
              LB;
              JB;
              Engle;
              sampleCorr;
              MASE;
              HitRate;
              stdErr];
