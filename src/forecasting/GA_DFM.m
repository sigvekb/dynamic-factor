function [measure] = ...
    GA_DFM(ga_vec)
%=========================
% Set up DFM environment
%=========================
g = 0;
max_iter = 12;
thresh = 1e-4;
selfLag = false;
restrictQ = false;
lags = 6;
dataFile = 'Oil_data.xlsx';
dataSheet = 'Oil_GA1';

horizons = (1);
outOfSampleMonths = 36;

% Set up block structure based on ga_vec
maxBlocks = 5;
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
% Set negative values to zero
blockStruct(blockStruct<0) = 0;
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
MASE = CalcMASE(DFMnorm(1:(end-oOSM),mainVar), forecasts_DFM(:,1)-mainActual);

[~,rmse] = CalcRMSFE(mainActual, forecasts_DFM(:,1), 1);

%measure = rmse;
measure = rmse+MASE;

if MASE < 1
    fileID = fopen('Runs.txt','a');
    fprintf(fileID,'\nMASE %2.3f, RMSE %2.3f\n',MASE, rmse);
    fprintf('\nMASE %2.3f, RMSE %2.3f\n',MASE, rmse);
    fclose(fileID);
    m = blockStruct(:);
    dlmwrite('Runs.txt',m', '-append', 'newline','pc');
end
