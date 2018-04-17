%=======================================
% DYNAMIC FACTOR MODEL W/ EM ALGORITHM
%=======================================
%*********************
% Input
%*********************
% directory - Set directory where you have data file
% dataFile - Name of excel file to read from
% sheet - Excel sheet name or number to read from
% block - The block structure of the data file
%         e.g. [3, 2] meaning first 3 columns make up first block,
%         and the 2 next columns make up second block
% maxIterations - Maximum number of iterations in EM algorithm
% lags - Number of lags, (for now restricted to be 1)

dir = 'C:\Users\Sigve Borgmo\OneDrive - NTNU\Indok\Master\dynamic-factor\Old Code';
dataFile = 'Dataset.xlsx';
dataSheet = 'Data';
blockFile = 'Blocks.xlsx';
blockSheet = 'B2';
outputFile = 'BlockSens';
maxIterations = 4000;
inflate = false;
threshold = -1;

writeRaw = true;
writeInput = true;
writeNormalized = true;
writeIMFIndex = true;
rng('default');

%*********************
% Preparation
%*********************
cd(dir);

% Data prep
[data, txt] = xlsread(dataFile, dataSheet, 'A3:BC273');
[inflateData, txt2] = xlsread(dataFile, 2);
preparedData = prepData(data, inflateData, inflate);
[T,n] = size(preparedData);

[blockStruct, txtB] = xlsread(blockFile, blockSheet, 'A1:I53');
[d,b]=size(blockStruct);
block = [];
selectedData = [];
rawData = [];
colName = txt(1,2:end);
colNewName = [];
for i=1:b
    block = [block, 0];
    for j=1:d
      commodity = blockStruct(j,i);
      if ~isnan(commodity)
          colNewName = [colNewName colName(commodity)];
          selectedData = [selectedData, preparedData(:,commodity)];
          rawData = [rawData, data(:,commodity)];
          block(i) = block(i)+1;
      end
    end
end

% Output prep
filename= strcat(outputFile,datestr(now,'yyyymmdd-HHMM'),'.xlsx');
factors = length(block)+1;
varNames = colNewName;
dates = txt(2:end-1,1)';
factorNames = ["Global", txtB];
% for i=1:length(block)
%     factorNames(i+1) = strcat("B", string(i));
% end
breakdown = ["Global", "Block", "Idio", "Total"];

%***********************
% Running the algorithm
%***********************
[normData,F_hat,F_pc,F_kal,num_iter, C, A, Q] = ...
    DynFactorDoz2(selectedData, factors, factors, 1, maxIterations, block);

%***********************
% Testing
%***********************
[var] = DynFactorTest(normData, F_hat, C);

%***********************
% Write to file
%***********************
xlswrite(filename,F_hat,"Factors",'B2');
xlswrite(filename,dates',"Factors",'A2');
xlswrite(filename,factorNames,"Factors",'B1');

xlswrite(filename,C,"Loadings",'B2');
xlswrite(filename,varNames',"Loadings",'A2');
xlswrite(filename,factorNames,"Loadings",'B1');

xlswrite(filename,var,"VarDecomp", 'B2');
xlswrite(filename,varNames',"VarDecomp",'A2');
xlswrite(filename,breakdown,"VarDecomp",'B1');

xlswrite(filename,A,"A",'B2');
xlswrite(filename,factorNames',"A",'A2');
xlswrite(filename,factorNames,"A",'B1');

xlswrite(filename,Q,"Q",'B2');
xlswrite(filename,factorNames,"Q",'B1');
xlswrite(filename,factorNames',"Q",'A2');

if writeRaw
    xlswrite(filename,rawData,"RawData",'B2');
    xlswrite(filename,varNames,"RawData",'B1');
    xlswrite(filename,dates',"RawData",'A2');
end

if writeInput
    xlswrite(filename,selectedData,"Input",'B2');
    xlswrite(filename,varNames,"Input",'B1');
    xlswrite(filename,dates',"Input",'A2');
end
if writeNormalized
    xlswrite(filename,normData,"NormInput",'B2');
    xlswrite(filename,varNames,"NormInput",'B1');
    xlswrite(filename,dates',"NormInput",'A2');
end

if writeIMFIndex
    title=["IMF PALLFNF", ""];
    xlswrite(filename,preparedData(:,54),"Index",'B2');
    xlswrite(filename,title,"Index",'B1');
    xlswrite(filename,dates',"Index",'A2');
end
