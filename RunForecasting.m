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
horizons = (1:3); %[1,3,6]
outOfSampleMonths = 36;

%===========
% DFM Input
%===========
blockFile = 'Blocks.xlsx';
blockSheet = 'Block2';

DFM = true;         % True: Run forecasting with DFM
globalFactors = 4;  % Number of global factors
maxIter = 400;      % Max number of iterations
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
blockStructure = blockData(2:end,2:end);

totalFactors = size(blockStructure,2)+globalFactors;
nanMatrix = CreateNaNMatrix(inputData);

[DFMData, newBlockStruct, DFMselection] = ...
            SelectData(inputData, blockStructure);


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


%=============
% Forecasting
%=============

% Run initial DFM
[~, factors, ~, A, C, Q, R] = ...
    DynamicFactorModel(DFMData, totalFactors, globalFactors, maxIter, threshold, ...
                       newBlockStruct, nanMatrix, lags, selfLag, restrictQ);

for t=1:outOfSampleMonths
    
    
end


