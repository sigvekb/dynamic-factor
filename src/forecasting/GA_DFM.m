function [measure] = ...
    GA_DFM(ga_vec)
%=========================
% Set up DFM environment
%=========================
g = 0;
max_iter = 10;
thresh = 1e-6;
selfLag = false;
restrictQ = false;
lags = 8;
dataFile = 'Salmon_data.xlsx';
dataSheet = 'GA_vars2';

horizons = (1);
outOfSampleMonths = 36;

% Set up block structure based on ga_vec
maxBlocks = 3;
maxVars = 3;

% Read data
[rawData, ~] = xlsread(dataFile, dataSheet, 'A1:FZ1000');
YoY = rawData(1,:);
LD = rawData(2,:);
rawData = rawData(3:end,:);

inputData = rawData;
inputData = LogDiff(true, inputData, YoY, LD);

% Salmon/oil must be the first variable in the dataset
blockStruct = [ones(1,maxBlocks); reshape(ga_vec, [maxVars, maxBlocks])];
oOSM = outOfSampleMonths;
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

VARlags = ones(1, size(blockStruct,2)+1)*lags;

[DFMData, newBlockStruct, ~] = ...
            SelectData(inputData, blockStruct);
        
% varNames = txt(1,2:end);
% varNames = varNames(DFMselection);

demean = bsxfun(@minus, DFMData, nanmean(DFMData));
DFMnorm = bsxfun(@rdivide, demean, nanstd(demean));
mainVar = 1;
mainActual = DFMnorm((end-oOSM+1):end,mainVar);

[forecasts_DFM, ~] = ForecastDFM(DFMData,horizons, oOSM, g, ...
                     max_iter, thresh, selfLag, restrictQ, newBlockStruct, VARlags, false);


% Calc wanted measure
naive_MAE = (1/(size(DFMnorm,1)-oOSM-1))*...
    nansum(abs(DFMnorm(2:(end-oOSM), mainVar)-DFMnorm(1:(end-oOSM-1), mainVar)));
forecastError = (1/oOSM)*nansum(abs(mainActual-forecasts_DFM(:,1)));
MASE = forecastError/naive_MAE;

[~,rmse] = CalcRMSFE(mainActual, forecasts_DFM(:,1), 1);

%measure = rmse;
measure = MASE;

if true
    fileID = fopen('Runs.txt','a');
    fprintf(fileID,'\n%2.3f\n',measure);
    fclose(fileID);
    m = blockStruct(:);
    dlmwrite('Runs.txt',m', '-append', 'newline','pc');
end
