%=============================================
% FORECASTING SCHEME W/ DFM + BENCHMARK MODELS
%=============================================
vars = 30;
maxBlocks = 5;
maxVars = 3;


lb = ones(1,maxVars*maxBlocks)*(-4);
ub = ones(1,maxVars*maxBlocks)*vars;
IntCon = (1:maxVars*maxBlocks);

solution = ga(@GA_DFM,maxVars*maxBlocks,[],[],[],[],lb,ub,[],IntCon);