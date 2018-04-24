%=============================================
% FORECASTING SCHEME W/ DFM + BENCHMARK MODELS
%=============================================

% dir = 'C:\Users\sigvekb\Master\dynamic-factor';
% dir = '\MATLAB Drive\Master\dynamic-factor';
dataFile = 'Dataset.xlsx';
dataSheet = 'Salmon2';
blockFile = 'Blocks.xlsx';
blockSheet = 'Block2';
% dataFile = 'SALMON.xlsm';
% dataSheet = 'Data1';
% blockFile = 'Block_WBC2.xlsx';
% blockSheet = 'Block2';
outputFile = 'DFM_Output';

globalFactors = 4;
maxIterations = 400;
threshold = 1e-6;
deflate = false;
logdiff = true;
selfLag = false; % Restrict factors to only load on own lags
restrictQ = false;