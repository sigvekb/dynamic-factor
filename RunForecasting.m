%=============================================
% FORECASTING SCHEME W/ DFM + BENCHMARK MODELS
%=============================================

%============
% Data input
%============
% dir = 'C:\Users\sigvekb\Master\dynamic-factor';
% dir = '\MATLAB Drive\Master\dynamic-factor';
dataFile = 'Salmon_data.xlsx';
dataSheet = 'Salmon2';
outputFile = 'ForecastingOutput';

%===================
% Forecasting input
%===================
horizons = [1,2,3];
outOfSampleMonths = 36;

%=========== 
% DFM Input
%===========
blockFile = 'Salmon_blocks.xlsx';
blockSheet = 'G5';

DFM = true;         % True: Run forecasting with DFM
globalFactors = 0;  % Number of global factors
maxIter = 1000;       % Max number of iterations
threshold = 1e-4;   % Convergence threshold for EM algorithm
deflate = false;    % True: Data is deflated according to US CPI
logdiff = true;     % True: Data is log differenced
selfLag = false;    % True: Restrict factors to only load on own lags
restrictQ = false;  % True: Q matrix is restricted to be diagonal

%======================
% Benchmark model input
%======================
modelFile = 'Benchmarks.xlsx';
modelSheet = 'Salmon';

ARIMA = true;
ARIMA_ar = 7;
ARIMA_ma = 2;

VAR = true;
VAR_lags = 4;

NOCHANGE = true;

%======================
% Data preparation
%======================
[rawData, txt] = xlsread(dataFile, dataSheet, 'A1:FZ1000');
YoY = rawData(1,:);
rawData = rawData(2:end,:);

inputData = rawData;
inputData = Deflate(deflate, inputData);
inputData = LogDiff(logdiff, inputData, YoY);

%==================
% DFM preparation
%==================
[blockData, blockTxt] = xlsread(blockFile, blockSheet, 'F1:AZ100');
lags = blockData(1,:);
blockStruct = blockData(2:end,2:end);

[DFMData, newBlockStruct, DFMselection] = ...
            SelectData(inputData, blockStruct);
        
demean = bsxfun(@minus, DFMData, nanmean(DFMData));
DFMnorm = bsxfun(@rdivide, demean, nanstd(demean));
mainVar = 1;
mainActual = DFMnorm((end-outOfSampleMonths+1):end,mainVar);

varNames = txt(1,2:end);
varNames = varNames(DFMselection);


%============================
% Benchark models preparation
%============================
[variables, models] = xlsread(modelFile, modelSheet, 'F1:AA100');

if VAR
    [VAR_data, VAR_struct, VAR_selection] = ...
                SelectData(inputData, variables(:,strcmp(models,'VAR')));
    demean = bsxfun(@minus, VAR_data, nanmean(VAR_data));
    VAR_norm = bsxfun(@rdivide, demean, nanstd(demean));
end

%=============
% Forecasting
%=============
if DFM
    [forecasts_DFM, varDecomp] = ...
        ForecastDFM(DFMData,horizons,outOfSampleMonths, globalFactors, ...
                    maxIter, threshold, selfLag, restrictQ, newBlockStruct, lags, true);
end
if VAR
    [forecasts_VAR] = ...
        ForecastVAR(VAR_norm, horizons, VAR_lags, outOfSampleMonths);
end
if ARIMA
    [forecasts_ARIMA] = ...
        ForecastARIMA(DFMnorm(:,mainVar), horizons, ARIMA_ar, ARIMA_ma, outOfSampleMonths);
end
if NOCHANGE
    [forecasts_NOCHANGE] = zeros(outOfSampleMonths, length(horizons));
    %ForecastARIMA(DFMnorm(:,mainVar), horizons, 0, 0, outOfSampleMonths);
end

%==============================
% Statistics to compare models
%==============================
if DFM
    mainForecast_DFM = permute(forecasts_DFM(:,mainVar,:),[1,3,2]);
    [Tf, H_len] = size(mainForecast_DFM);
    benchmarks = zeros(Tf,H_len, 3);
    if VAR
        benchmarks(:,:,1) = permute(forecasts_VAR(:,mainVar,:),[1,3,2]);
    end
    if ARIMA
        benchmarks(:,:,2) = permute(forecasts_ARIMA(:,mainVar,:),[1,3,2]);
    end
    if NOCHANGE
        benchmarks(:,:,3) = forecasts_NOCHANGE;
    end
    [statistics] = ForecastStatistics(mainActual, mainForecast_DFM, ...
                                        benchmarks, horizons);
end

%====================
% Creating dashboard
%====================
M=3;
N=H_len;

p = panel();
p.pack(M, N)
p.de.margin = 8;
p.margin = [7 7 7 7];

clf;
error = mainActual - mainForecast_DFM(:,1);
maxErr = max(abs(error));
for h=1:N
    %==========
    % PANEL 1
    %==========
    p(1,h).select();
    plot(mainActual,'k','LineWidth',1.5);
    hold on;
    plot(mainForecast_DFM(:,h),'r-','LineWidth',1.5);
    if VAR
        plot(benchmarks(:,h,1),'b-','LineWidth',1.0);
    end
    if ARIMA
        plot(benchmarks(:,h,2),'g-','LineWidth',1.0);
    end
    if NOCHANGE
        plot(benchmarks(:,h,3),'y-','LineWidth',1.0);
    end

    if h==1
        legend('Actual','DFM','VAR','ARIMA','NOCHANGE','Location','NorthWest');
    end
    
    title(strcat('Forecast Horizon, h = ',num2str(horizons(h))));
    axis tight;
    set(gca,'FontSize',10,'FontName','Times');
    grid on;
    ax = gca;
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.2;
    ax.FontSize = 10;
    ax.LineWidth = 0.5;
    
    %==========
    % PANEL 2
    %==========
    p(2,h).select();
    
    r_VAR = statistics(3,h);
    r_ARIMA = statistics(4,h);
    r_naive = statistics(5,h);

    VARstars = getStars(statistics(6,h));
    ARIMAstars = getStars(statistics(7,h));
    naivestars = getStars(statistics(8,h));
    LBstars = getStars(statistics(9,h));
    JBstars = getStars(statistics(10,h));
    ENstars = getStars(statistics(11,h));
    
    rmse_txt_VAR = strcat('rRMSE_{VAR} = ', num2str(r_VAR, '%1.3f'), VARstars);
    rmse_txt_ARIMA = strcat('rRMSE_{ARIMA} = ',num2str(r_ARIMA,'%1.3f'), ARIMAstars);
    rmse_txt_NOCHANGE = strcat('rRMSE_{NC} = ',num2str(r_naive,'%1.3f'), naivestars);
    lb_txt = strcat('Ljung-Box:', num2str(statistics(9,h),'%1.2f'), LBstars);
    jb_txt = strcat('Jarque Bera:', num2str(statistics(10,h),'%1.2f'), JBstars);
    en_txt = strcat('Engle Test:', num2str(statistics(11,h),'%1.2f'), ENstars);
    corr_txt = strcat('Sample Corr:',num2str(statistics(12,h),'%1.2f'));
    relCorr_VAR = strcat('rCorr_{VAR} = ', num2str(statistics(13,h), '%1.3f'));
    relCorr_ARIMA = strcat('rCorr_{ARIMA} = ', num2str(statistics(14,h), '%1.3f'));
    relCorr_NOCHANGE = strcat('rCorr_{NC} = ', num2str(statistics(15,h), '%1.3f'));
    
    % Text box
    x = -0.3+0.33*h;
    dim1 = [x 0.2 0.2 0.42]; %location of box on the panel
    y = x + 0.1666;
    dim2 = [y 0.2 0.2 0.42];
    annotation(...
        'textbox',dim1,...
        'String', {'Benchmark statistics:',' ', ...
                   rmse_txt_VAR,rmse_txt_ARIMA,rmse_txt_NOCHANGE,...
                   relCorr_VAR,relCorr_ARIMA,relCorr_NOCHANGE},...
        'FitBoxToText','on',...
        'FontSize',10,...
        'FontName','Times');
    annotation(...
        'textbox',dim2,...
        'String', {'Forecast error statistics:',' ', ...
                   lb_txt,jb_txt,en_txt,corr_txt},...
        'FitBoxToText','on',...
        'FontSize',10,...
        'FontName','Times');

    hold on; 

    %==========
    % PANEL 3
    %==========
    p(3,h).select();
    lower_ci = mainForecast_DFM(:,h)-1.96*statistics(13,h);
    upper_ci = mainForecast_DFM(:,h)+1.96*statistics(13,h); 

    h1 = plot(lower_ci','k--','LineWidth',1);
    hold on;
    h2 = plot(upper_ci','k--','LineWidth',1);

    plot(mainActual,'k','LineWidth',1.5);
    plot(mainForecast_DFM(:,h),'r','LineWidth',3);

    z = (1:outOfSampleMonths)';
    fill( [z' fliplr(z')],  [lower_ci' fliplr(upper_ci')], 'r', 'EdgeColor','none');
    alpha(0.1);
    
    if h==1
        legend('95% CI','Location','NorthWest');
    end

    set(gca,'FontSize',10,'FontName','Times');
    axis tight;
    grid on;
    ax = gca;
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.2;
    ax.FontSize = 10;
    ax.LineWidth = 0.5;
end