%=======================================
% DYNAMIC FACTOR MODEL W/ EM ALGORITHM
%=======================================
%*********************
% Input
%*********************
% dir - Set directory where you have data file
% dataFile - Name of excel file to read from
% dataSheet - Excel sheet name or number to read from
% blockFile - Name of Excel file containing block structure
% blockSheet - Excel sheet name or number to read block structure from
% maxIterations - Maximum number of iterations in EM algorithm
% threshold - The threshold value for convergence of the EM algorithm.
%               Should be between 1e-4 and 1e-7

% dir = 'C:\Users\sigvekb\Master\dynamic-factor';
% dir = '\MATLAB Drive\Master\dynamic-factor';
dataFile = 'Dataset.xlsx';
dataSheet = 'Comm2';
blockFile = 'Blocks.xlsx';
blockSheet = 'Comm2';
% dataFile = 'SALMON.xlsm';
% dataSheet = 'Data1';
% blockFile = 'Block_WBC2.xlsx';
% blockSheet = 'Block2';
outputFile = 'DFM_Output';

globalFactors = 1;
maxIterations = 100;
threshold = 1e-7;
deflate = false;
logdiff = true;
selfLag = true; % Restrict factors to only load on own lags
restrictQ = true;

writeRaw = false;
writeInput = false;
writeNormalized = false;

%cd(dir);
%==================
% DFM preparation
%==================
[rawData, txt] = xlsread(dataFile, dataSheet, 'A1:FZ1000');
YoY = rawData(1,:);
rawData = rawData(2:end,:);

inputData = rawData;
inputData = Deflate(deflate, inputData);
inputData = LogDiff(logdiff, inputData, YoY);

[blockData, blockTxt] = xlsread(blockFile, blockSheet, 'F1:AZ100');
lags = blockData(1,:);
blockStruct = blockData(2:end,2:end);

[DFMData, newBlockStruct, selection] = ...
            SelectData(inputData, blockStruct);

%=======================
% Running the algorithm
%=======================
[normData, F_hat, iter, C, A, Q] = ...
    DynamicFactorModel(DFMData, totalFactors, globalFactors, maxIterations, threshold, ...
                       newBlockStruct, nanMatrix, lags, selfLag, restrictQ, ...
                       [],[],[],[]);

fprintf('Finished in %d iterations', iter-1);

%***********************
% Variance Decomposition
%***********************
[varDecomp] = VarianceDecomposition(normData, F_hat, C, globalFactors);

%***********************
% Write to file
%***********************
% Output prep
dates = txt(2:end-1,1)';
varNames = txt(1,2:end);
varNames = varNames(selection);
rawData = rawData(:, selection);
factorNames = blockTxt(1,2:end);
factorNames = [repelem(("Global"),globalFactors) factorNames];
maxlags = max(lags);

outputFile = strcat(outputFile,datestr(now,'mmdd-HHMM'),'.xlsx');
breakdown = [factorNames 'Global' 'Block' 'Idio'];

xlswrite(outputFile,F_hat,'Factors','B2');
xlswrite(outputFile,dates','Factors','A2');
xlswrite(outputFile,factorNames,'Factors','B1');

xlswrite(outputFile,C,'Loadings','B2');
xlswrite(outputFile,varNames','Loadings','A2');
xlswrite(outputFile,factorNames,'Loadings','B1');

xlswrite(outputFile,varDecomp,'VarDecomp', 'B2');
xlswrite(outputFile,varNames','VarDecomp','A2');
xlswrite(outputFile,breakdown,'VarDecomp','B1');

AfactorNames = repmat(factorNames, 1, maxlags);
xlswrite(outputFile,A,'A','B2');
xlswrite(outputFile,factorNames','A','A2');
xlswrite(outputFile,factorNames,'A','B1');

xlswrite(outputFile,Q,'Q','B2');
xlswrite(outputFile,factorNames,'Q','B1');
xlswrite(outputFile,factorNames','Q','A2');

if writeRaw
    xlswrite(outputFile,rawData,'RawData','B2');
    xlswrite(outputFile,varNames,'RawData','B1');
    xlswrite(outputFile,dates','RawData','A2');
end

if writeInput
    xlswrite(outputFile,preparedData,'Input','B2');
    xlswrite(outputFile,varNames,'Input','B1');
    xlswrite(outputFile,dates','Input','A2');
end

if writeNormalized
    xlswrite(outputFile,normData,'NormInput','B2');
    xlswrite(outputFile,varNames,'NormInput','B1');
    xlswrite(outputFile,dates','NormInput','A2');
end