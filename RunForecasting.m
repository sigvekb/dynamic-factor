%=============================================
% FORECASTING SCHEME W/ DFM + BENCHMARK MODELS
%=============================================

%============
% Data input
%============
% dir = 'C:\Users\sigvekb\Master\dynamic-factor';
% dir = '\MATLAB Drive\Master\dynamic-factor';
dataFile = 'Dataset.xlsx';
dataSheet = 'Salmon2';
outputFile = 'ForecastingOutput';

%===================
% Forecasting input
%===================
horizons = [1,3,6];
outOfSampleMonths = 36;

%===========
% DFM Input
%===========
blockFile = 'Blocks1.xlsx';
blockSheet = 'B6';

DFM = true;         % True: Run forecasting with DFM
globalFactors = 0;  % Number of global factors
maxIter = 100;       % Max number of iterations
threshold = 1e-6;   % Convergence threshold for EM algorithm
deflate = false;    % True: Data is deflated according to US CPI
logdiff = true;     % True: Data is log differenced
selfLag = false;    % True: Restrict factors to only load on own lags
restrictQ = false;  % True: Q matrix is restricted to be diagonal

%======================
% Benchmark model input
%======================
modelFile = 'Benchmarks.xlsx';
modelSheet = 'VAR';

ARIMA = false;
ARIMA_ar = 3;
ARIMA_ma = 2;

VAR = true;
VAR_lags = 4;

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

varNames = txt(1,2:end);
varNames = varNames(DFMselection);


%============================
% Benchark models preparation
%============================
[variables, models] = xlsread(modelFile, modelSheet, 'F1:AA100');

% Create datasets
if VAR
    [VAR_data, VAR_struct, VAR_selection] = ...
                SelectData(inputData, variables(:,strcmp(models,'VAR')));
    demean = bsxfun(@minus, VAR_data, nanmean(VAR_data));
    VAR_norm = bsxfun(@rdivide, demean, nanstd(demean));
end
if ARIMA
    [ARIMA_data, ARIMA_struct, ARIMA_selection] = ...
                SelectData(inputData, variables(:,strcmp(models,'ARIMA')));
    demean = bsxfun(@minus, ARIMA_data, nanmean(ARIMA_data));
    ARIMA_norm = bsxfun(@rdivide, demean, nanstd(demean));
end


%========================
% Forecasting preparation
%========================
h_periods = length(horizons);
maxHorizon = max(horizons);

[T,n] = size(DFMData);
DFM_forecasts = zeros(outOfSampleMonths,n,h_periods);

nVAR = size(VAR_norm,2);
VAR_forecasts = zeros(outOfSampleMonths,nVAR,h_periods);
VARModel = varm(nVAR, VAR_lags);

%==============
% Forecasting
%==============
for t=1:outOfSampleMonths
    fprintf('\nForecasting month: %d of %d\n', t, outOfSampleMonths);
    removeMonths = outOfSampleMonths-t+1;
    
    % Forecast DFM
    if DFM
        
        forecastData = [DFMData(1:(end-removeMonths),:); NaN([maxHorizon,n])];
        [~, factors, ~, A, C, Q, R,initV] = ...
                DynamicFactorModel(forecastData, globalFactors, maxIter, ...
                    threshold, selfLag, restrictQ, newBlockStruct, lags, []);

        for h=1:h_periods
            horizon=horizons(h);
            if removeMonths>=horizon
                f_index = size(factors,1)-maxHorizon+horizon;
                y_f = C*factors(f_index,:)';
                DFM_forecasts(horizon+t-1,:,h) = y_f';
            end
        end
    end

    % Forecast VAR
    if VAR
        forecastData_VAR = VAR_norm(1:(end-removeMonths),:);
        
        VARFit = estimate(VARModel,forecastData_VAR,'display','off');
        VARForecasts = forecast(VARFit,maxHorizon,forecastData_VAR);
        
        for h=1:h_periods
            horizon=horizons(h);
            if removeMonths>=horizon
                VAR_forecasts(horizon+t-1,:,h) = VARForecasts(horizon, :);
            end
        end
    end
end


if DFM
    % Variance Decomposition of DFM
    demean = bsxfun(@minus, DFMData, nanmean(DFMData));
    DFMnorm = bsxfun(@rdivide, demean, nanstd(demean));  
    [varDecomp] = VarianceDecomposition(DFMnorm, factors(1:(end-maxHorizon+1),:), C, globalFactors);
end

%==================================
% Calculate forecast error measures
%==================================
if DFM
    DFM_forecasts(DFM_forecasts == 0) = NaN;
    RMSE_DFM = zeros(n,h_periods);
    RMSFE_DFM = zeros(n,h_periods);
    actual_DFM = DFMnorm(T-outOfSampleMonths+1:end,:);
end
if VAR
    VAR_forecasts(VAR_forecasts == 0) = NaN;
    RMSE_VAR = zeros(nVAR,h_periods);
    RMSFE_VAR = zeros(nVAR,h_periods);
    RMSE_relative = zeros(h_periods,1);
    DM = zeros(h_periods,1); % Retrieve P-value from Diebold-Mariano
    actual_VAR = VAR_norm(T-outOfSampleMonths+1:end,:);
end


for i=1:h_periods
    h = horizons(i);
    if DFM
        forecast_DFM = DFM_forecasts(:,:,i);
        [rmsfe, rmse] = RMSFE(actual_DFM, forecast_DFM, h);
        RMSE_DFM(:,i) = rmse';
        RMSFE_DFM(:,i) = rmsfe';
    end
    if VAR
        forecast_VAR = VAR_forecasts(:,:,i);
        [rmsfe, rmse] = RMSFE(actual_VAR, forecast_VAR, h);
        RMSE_VAR(:,i) = rmse';
        RMSFE_VAR(:,i) = rmsfe';
    end
end

%==============================
% Statistics to compare models
%==============================
% Calculate relative RMSE
salmonIndex = 1;
if DFM && VAR
    for i=1:h_periods
        h = horizons(i);
        RMSE_relative(i) = RMSFE_DFM(salmonIndex,i) / RMSFE_VAR(salmonIndex,i);
        forecast_DFM = DFM_forecasts(h:end,1,i);
        forecast_VAR = VAR_forecasts(h:end,1,i);
        error_DFM = actual_DFM(h:end,salmonIndex)-forecast_DFM;
        error_VAR = actual_VAR(h:end,salmonIndex)-forecast_VAR;
        [~,p] = DieboldMariano(error_DFM, error_VAR);
        DM(i) = p;
    end
end

%=======
% Plots
%=======

salmonForecast = permute(DFM_forecasts(:,1,:),[1,3,2]);
salmonForecast_VAR = permute(VAR_forecasts(:,1,:),[1,3,2]);
salmonActual = DFMnorm((T-outOfSampleMonths+1):end,1);

%====PLOT===
M=3;
N=h_periods;

p = panel();
p.pack(M, N)

% panel marigns
p.de.margin = 8;
p.margin = [6 6 6 6];

% Antall kolonner må inn manuelt (M) !
clf;
error = salmonActual - salmonForecast(:,1);
maxErr = max(abs(error));
for i = 1:N
        p(1,i).select();
        plot(salmonActual,'k','LineWidth',2);
        hold on;
        plot(salmonForecast(:,i),'r--','LineWidth',1);
        plot(salmonForecast_VAR(:,i),'b--','LineWidth',1)
        legend('PSALM','DFM Forecast','Location','NorthWest');
        step = num2str(i);
        txt = strcat('Forecast horizon h = ',step);
        title(txt);
        axis tight;
        
        set(gca,'FontSize',10,'FontName','Times');
        grid on;
        ax = gca;
        ax.GridLineStyle = ':';
        ax.GridAlpha = 0.2;
        ax.FontSize = 10;
        ax.LineWidth = 0.5;
        %==========================
        %== FORECAST ERROR ABS=====
        %=========================
        p(2,i).select();
        error = salmonActual - salmonForecast(:,i);
        bar(1:outOfSampleMonths,abs(error),'r');
        
        set(gca,'FontSize',10,'FontName','Times');
        grid on;
        ax = gca;
        ax.GridLineStyle = ':';
        ax.GridAlpha = 0.2;
        ax.FontSize = 10;
        ax.LineWidth = 0.5;
        
        axis tight;
        ylim([0,maxErr+0.2*maxErr]);
        hold on; 
        
        %=============================
        %=== FORECASTR ERROR SCATTER==
        %=============================
        p(3,i).select();
        scatter(1:outOfSampleMonths,error,'filled','r');
        
        set(gca,'FontSize',10,'FontName','Times');
        grid on;
        ax = gca;
        ax.GridLineStyle = ':';
        ax.GridAlpha = 0.2;
        ax.FontSize = 10;
        ax.LineWidth = 0.5;
        
        r = RMSE_relative(i);
        rmse_num = num2str(r, '%1.3f');
        pVal1 = DM(i,1);
        if pVal1 < 0.01
            DMstars = '***';
        elseif pVal1 < 0.05
            DMstars = '**';
        elseif pVal1 < 0.10
            DMstars = '*';
        else
            DMstars = '';
        end
        rmse_txt = strcat('RMSE DFM/VAR = ', rmse_num, DMstars);
%         r = RMSE_relative_naive(i);
%         rmse_num = num2str(r, '%1.3f');
%         pVal2 = DM(i,2);
%         if pVal2 < 0.01
%             DMstars = '***';
%         elseif pVal2 < 0.05
%             DMstars = '**';
%         elseif pVal2 < 0.10
%             DMstars = '*';
%         else
%             DMstars = '';
%         end
%         rmse_naive_text = strcat('RMSE DFM/Naive = ', rmse_num, DMstars);
        
        x = 0.03+(1/h_periods)*(i-1);
        dim = [x 0.2 0.2 0.44];
        annotation(...
            'textbox',dim,...
            'String', {rmse_txt,strcat('p=',num2str(pVal1,'%1.3f'))},...
            'FitBoxToText','on',...
            'FontSize',8,...
            'FontName','Times');
        
        ylim([-(maxErr+0.2),maxErr+0.2]);
        xlim([0,outOfSampleMonths]);
        
        hline = refline(0,0);
        hline.Color = 'k';
        hline.LineWidth = 0.01;
        hline.LineStyle = ':';
        hold on;    
end
