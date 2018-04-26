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
horizons = (1:2); %[1,3,6]
outOfSampleMonths = 6;

%===========
% DFM Input
%===========
blockFile = 'Blocks.xlsx';
blockSheet = 'Block2';

DFM = true;         % True: Run forecasting with DFM
globalFactors = 1;  % Number of global factors
maxIter = 25;       % Max number of iterations
threshold = 1e-6;   % Convergence threshold for EM algorithm
deflate = false;    % True: Data is deflated according to US CPI
logdiff = true;     % True: Data is log differenced
selfLag = false;    % True: Restrict factors to only load on own lags
restrictQ = false;  % True: Q matrix is restricted to be diagonal

%======================
% Benchmark model input
%======================
modelFile = 'Forecasting.xlsx';
modelSheet = 'Model';

ARIMA = true;
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

totalFactors = size(blockStruct,2)+globalFactors;
nanMatrix = CreateNaNMatrix(inputData);

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
    VAR_index = find(strcmp(models,'VAR'));
    [VAR_data, VAR_struct, VAR_selection] = ...
                SelectData(inputData, variables(:,VAR_index));
end
if ARIMA
    ARIMA_index = find(strcmp(models,'ARIMA'));
    [ARIMA_data, ARIMA_struct, ARIMA_selection] = ...
                SelectData(inputData, variables(:,ARIMA_index));
end


%========================
% Forecasting preparation
%========================
% Introduce necessary variables/collectors, such that all relevant data
% from the forecasting scheme is saved
h_periods = length(horizons);
maxHorizon = max(horizons);

[T,n] = size(DFMData);
DFM_forecasts = zeros(outOfSampleMonths,n,h_periods);

%=============
% Forecasting
%=============
% Run initial DFM
% [DFMnorm, originalFactors, ~, A, C, Q, R,initV] = ...
%     DynamicFactorModel(DFMData, globalFactors, maxIter, ...
%                 threshold, selfLag, restrictQ, blockStruct, lags, []);
            
demean = bsxfun(@minus, DFMData, nanmean(DFMData));
DFMnorm = bsxfun(@rdivide, demean, nanstd(demean));
originalC = C;    

for t=1:outOfSampleMonths
    fprintf("Forecast: %d", t);
    % Forecast DFM
    removeMonths = outOfSampleMonths-t+1;
    forecastData = [DFMData(1:(end-removeMonths),:); NaN([maxHorizon,n])];
    [~, factors, ~, A, C, Q, R,initV] = ...
            DynamicFactorModel(forecastData, globalFactors, maxIter, ...
                threshold, selfLag, restrictQ, newBlockStruct, lags, []);

    for h=1:h_periods
        horizon=horizons(h);
        if removeMonths>=horizon
            f_index = size(factors,1)-maxHorizon+horizon;
            forecast = C*factors(f_index,:)';
            DFM_forecasts(horizon+t-1,:,h) = forecast';
        end
    end
end

[varDecomp] = VarianceDecomposition(DFMnorm, factors(1:(end-maxHorizon+1),:), C, globalFactors);

% Set non-forecasted values to NaN
DFM_forecasts(DFM_forecasts == 0) = NaN;

RMSE = zeros(n,h_periods);
RMSFE = zeros(n,h_periods); % RMSE that accounts for horizon length 
actualValues = DFMnorm(T-outOfSampleMonths+1:end,:);


for i=1:h_periods
    horizon = horizons(i);
    forecastValues = DFM_forecasts(:,:,i);
    mse = mean((actualValues-forecastValues).^2,1,'omitnan');
    rmse = sqrt(mse);
    RMSE(:,i) = rmse';
    rmsfe = sqrt(mse/horizon);
    RMSFE(:,i) = rmsfe';
end

salmonForecast = permute(DFM_forecasts(:,8,:),[1,3,2]);
salmonActual = DFMnorm((end-outOfSampleMonths+1):end,8);
