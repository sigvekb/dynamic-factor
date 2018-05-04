function [inverse] = CholeskyInversion(matrix)
try
    choleskyDecomp = chol(matrix);
catch ME
    warning('Matrix is not PD, nearestSPD is used instead');
    mHat = nearestSPD(matrix);
    choleskyDecomp = chol(mHat);
end
I = eye(size(matrix,1));
inverse = I / choleskyDecomp / choleskyDecomp';