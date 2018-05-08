%=============================================
% FORECASTING SCHEME W/ DFM + BENCHMARK MODELS
%=============================================
addpath(genpath('dynamic-factor/src'))

vars = 23;
maxVars = 3;
maxBlocks = 3;

lb = zeros(1,maxVars*maxBlocks);
ub = ones(1,maxVars*maxBlocks)*vars;
IntCon = (1:maxVars*maxBlocks);

solution = ga(@GA_DFM,maxVars*maxBlocks,[],[],[],[],lb,ub,[],IntCon);