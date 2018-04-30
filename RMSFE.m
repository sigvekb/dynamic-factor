function [rmsfe, rmse] = RMSFE(y, y_f, h)
    mse = mean((y-y_f).^2,1,'omitnan');
    rmse = sqrt(mse);
    rmsfe = sqrt(mse/h);

