%=============================================
% FORECASTING SCHEME W/ DFM + BENCHMARK MODELS
%=============================================

%============
% Data input
%============
% dataFile = 'Salmon_data.xlsx';
% dataSheet = 'Salmon2';
dataFile = 'Oil_data.xlsx';
dataSheet = 'Raw2006';

%===================
% Forecasting input
%===================
horizons = [1,2,3];
outOfSampleMonths = 36;

%===========
% DFM Input
%===========
% blockFile = 'Salmon_blocks.xlsx';
% blockSheet = 'G5';
blockFile = 'Oil_blocks.xlsx';
blockSheet = 'GA3';

DFM = true;         % True: Run forecasting with DFM
globalFactors = 0;  % Number of global factors
maxIter = 12;       % Max number of iterations
threshold = 1e-5;   % Convergence threshold for EM algorithm
deflate = false;    % True: Data is deflated according to US CPI
logdiff = true;     % True: Data is log differenced
selfLag = false;    % True: Restrict factors to only load on own lags
restrictQ = false;  % True: Q matrix is restricted to be diagonal

%======================
% Benchmark model input
%======================
modelFile = 'Benchmarks.xlsx';
modelSheet = 'Oil';

ARIMA = true;
ARIMA_ar = 5; % Salmon ARIMA(7,0,2) Oil ARIMA(5,0,2)
ARIMA_ma = 2;

VAR = true;
VAR_lags = 3; % Salmon VAR(4), Oil VAR(1)

NOCHANGE = true;

%======================
% Data preparation
%======================
[rawData, txt] = xlsread(dataFile, dataSheet, 'A1:FZ1000');
YoY = rawData(1,:);
LD = rawData(2,:);
rawData = rawData(3:end,:);

inputData = rawData;
inputData = Deflate(deflate, inputData);
inputData = LogDiff(logdiff, inputData, YoY, LD);

%==================
% DFM preparation
%==================
[blockData, blockTxt] = xlsread(blockFile, blockSheet, 'F1:AZ200');
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


%=============================
% Benchark models preparation
%=============================
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
    [forecasts_DFM, varDecomp, C] = ...
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
    %models(:,:,1) = mainForecast_DFM;
    models = zeros(Tf,H_len, 3);
    if VAR
        models(:,:,3) = permute(forecasts_VAR(:,mainVar,:),[1,3,2]);
    end
    if ARIMA
        models(:,:,2) = permute(forecasts_ARIMA(:,mainVar,:),[1,3,2]);
    end
    if NOCHANGE
        models(:,:,1) = forecasts_NOCHANGE;
    end
    [statistics] = ForecastStatistics(mainActual, mainForecast_DFM, ...
                         models, horizons, DFMnorm,outOfSampleMonths);
end
error_DFM = bsxfun(@minus, mainForecast_DFM, mainActual);

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
        plot(models(:,h,1),'b-','LineWidth',1.0);
    end
    if ARIMA
        plot(models(:,h,2),'b--','LineWidth',1.0);
    end
    if NOCHANGE
        plot(models(:,h,3),'b:','LineWidth',1.0);
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
    
    r_VAR = statistics(5,h);
    r_ARIMA = statistics(6,h);
    r_naive = statistics(7,h);

    VARstars = getStars(statistics(8,h));
    ARIMAstars = getStars(statistics(9,h));
    naivestars = getStars(statistics(10,h));
    LBstars = getStars(statistics(11,h));
    JBstars = getStars(statistics(15,h));
    ENstars = getStars(statistics(19,h));
    
    rmse_txt_VAR = strcat('rRMSE_{VAR} = ', num2str(r_VAR, '%1.3f'), VARstars);
    rmse_txt_ARIMA = strcat('rRMSE_{ARIMA} = ',num2str(r_ARIMA,'%1.3f'), ARIMAstars);
    rmse_txt_NOCHANGE = strcat('rRMSE_{NC} = ',num2str(r_naive,'%1.3f'), naivestars);
    lb_txt = strcat('Ljung-Box:', num2str(statistics(11,h),'%1.2f'), LBstars);
    jb_txt = strcat('Jarque Bera:', num2str(statistics(15,h),'%1.2f'), JBstars);
    en_txt = strcat('Engle Test:', num2str(statistics(19,h),'%1.2f'), ENstars);
    corr_txt = strcat('Sample Corr:',num2str(statistics(23,h),'%1.2f'));
    relCorr_NOCHANGE = strcat('Corr_{NC} = ', num2str(statistics(28,h), '%1.3f'));
    relCorr_ARIMA = strcat('Corr_{ARIMA} = ', num2str(statistics(29,h), '%1.3f'));
    relCorr_VAR = strcat('Corr_{VAR} = ', num2str(statistics(30,h), '%1.3f'));
    MASE_txt = strcat('MASE = ', num2str(statistics(27,h), '%1.3f'));
    HR_txt = strcat('HitRate = ', num2str(statistics(31,h), '%1.3f'));
    iter_txt = strcat('Iterations = ', num2str(maxIter, '%3d'));
    maxlag_txt = strcat('Lag = ', num2str(mode(lags), '%2d'));
    
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
                   lb_txt,jb_txt,en_txt,corr_txt, MASE_txt,HR_txt,iter_txt,maxlag_txt},...
        'FitBoxToText','on',...
        'FontSize',10,...
        'FontName','Times');

    hold on; 

    %==========
    % PANEL 3
    %==========
    p(3,h).select();
    lower_ci = mainForecast_DFM(:,h)-1.96*statistics(32,h);
    upper_ci = mainForecast_DFM(:,h)+1.96*statistics(32,h); 

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