function [HR] = CalcHR(y,y_hat)
% Calculate Hit Rate for forecast versus actual   
yTimesyHat = y.*y_hat;
s=sign(yTimesyHat);
numPos=sum(s(:)==1);
HR = numPos/size(y_hat,1);