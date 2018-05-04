function [rmse] = ...
    GA_DFM(ga_vec)
%=========================
% Set up DFM environment
%=========================
g = 0;
max_iter = 10;
thresh = 1e-6;
selfLag = false;
restrictQ = false;
lags = 4;
dataFile = '..\..\data\Salmon_data.xlsx';
dataSheet = 'GA_vars';

horizons = [1];
outOfSampleMonths = 24;

% Read data
[rawData, ~] = xlsread(dataFile, dataSheet, 'A1:FZ1000');
YoY = rawData(1,:);
rawData = rawData(2:end,:);

inputData = rawData;
inputData = LogDiff(true, inputData, YoY);

% Set up block structure based on ga_vec
maxBlocks = 2;
maxVars = 2;

% Salmon/oil must be the first variable in the dataset
blockStruct = [ones(1,maxBlocks); reshape(ga_vec, [maxVars, maxBlocks])];

% Set duplicate values to zero
remove = [];
for i=1:maxBlocks
    [~, ind] = unique(blockStruct(:,i), 'rows');
    dup_ind = setdiff(1:(maxVars+1), ind);
	% duplicate values
    blockStruct(dup_ind, i) = 0;
    if sum(blockStruct(:,i) == 0) == maxVars
        remove = [remove i];
    end
end
blockStruct(:,remove) = [];
blockStruct(2:end,:) = sort(blockStruct(2:end,:), 'descend');
blockStruct(blockStruct == 0) = NaN;

VARlags = ones(1, size(blockStruct,1)+1)*lags;

[DFMData, newBlockStruct, ~] = ...
            SelectData(inputData, blockStruct);

demean = bsxfun(@minus, DFMData, nanmean(DFMData));
DFMnorm = bsxfun(@rdivide, demean, nanstd(demean));
mainVar = 1;
mainActual = DFMnorm((end-outOfSampleMonths+1):end,mainVar);

[forecasts_DFM, ~] = ForecastDFM(DFMData,horizons, outOfSampleMonths, g, ...
                     max_iter, thresh, selfLag, restrictQ, newBlockStruct, VARlags, false);


[~,rmse] = CalcRMSFE(mainActual, forecasts_DFM(:,1,1), 1);