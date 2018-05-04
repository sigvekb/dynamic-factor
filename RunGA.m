addpath(genpath('dynamic-factor/src'))

maxVars = 2;
maxBlocks = 2;
lb = [0 0 0 0];
ub = ones(1,4)*10;
% %ub = [19 25 25 25,...
%       25 25 25 25,...
%       25 25 25 25,...
%       25 25 25 25];
IntCon = (1:4);

solution = ga(@GA_DFM,maxVars*maxBlocks,[],[],[],[],lb,ub,[],IntCon);

a = 1;