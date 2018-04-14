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
% lags - Number of lags, (for now restricted to be 1)
% threshold - The threshold value for convergence of the EM algorithm.
%               Should be between 1e-4 and 1e-7

dir = '/Users/sigvekb/Master/dynamic-factor';
dataFile = 'WorldBankCommodities.xlsx';
dataSheet = 'Data';
blockFile = 'Block_WBC.xlsx';
blockSheet = 'B1';
outputFile = 'DFM_Output';
maxIterations = 1000;
threshold = 1e-5;

deflate = false;
logdiff = true;

writeRaw = true;
writeInput = false;
writeNormalized = false;
writeIMFIndex = false;

%*********************
% Preparation
%*********************
cd(dir);

% Data preparation
[data, txt]                  = xlsread(dataFile, dataSheet, 'A1:BE460');
[blockStructure, blockNames] = xlsread(blockFile, blockSheet, 'A1:I53');

[preparedData, nanMatrix, blockCount, selection] = ... 
            PrepareData(data, deflate, logdiff, blockStructure);

r = length(blockCount)+1;

%***********************
% Running the algorithm
%***********************
[normData, F_hat, iter, C, A, Q] = ...
    DynamicFactorModel(preparedData, r, r, 1, maxIterations, ...
                       threshold, blockCount, nanMatrix);

fprintf('Finished in %d iterations', iter-1);

%***********************
% Testing
%***********************
[var] = DynFactorTest(normData, F_hat, C);

%***********************
% Write to file
%***********************
% Output prep
dates = txt(2:end-1,1)';
varNames = txt(1,2:end);
varNames = varNames(selection);
rawData = data(:, selection);

outputFile = strcat(outputFile,datestr(now,'mmdd-HHMM'),'.xlsx');
factorNames = ['Global', blockNames];
breakdown = ["Global", "Block", "Idio", "Total"];

xlswrite(outputFile,F_hat,'Factors','B2');
xlswrite(outputFile,dates','Factors','A2');
xlswrite(outputFile,factorNames,'Factors','B1');

xlswrite(outputFile,C,'Loadings','B2');
xlswrite(outputFile,varNames','Loadings','A2');
xlswrite(outputFile,factorNames,'Loadings','B1');

xlswrite(outputFile,var,'VarDecomp', 'B2');
xlswrite(outputFile,varNames','VarDecomp','A2');
xlswrite(outputFile,breakdown,'VarDecomp','B1');

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

if writeIMFIndex
    title=['IMF PALLFNF', ''];
    xlswrite(outputFile,preparedData(:,54),'Index','B2');
    xlswrite(outputFile,title,'Index','B1');
    xlswrite(outputFile,dates','Index','A2');
end
