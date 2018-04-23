function [inverse] = CholeskyInversion(matrix)
choleskyDecomp = chol(matrix);
I = eye(size(matrix,1));
inverse = I / choleskyDecomp / choleskyDecomp';