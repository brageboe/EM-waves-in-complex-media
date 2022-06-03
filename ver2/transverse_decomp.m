function [matrix_tt, matrix_t, matrix_z, matrix_zz] = transverse_decomp(matrix)
% decomposition as in lecture notes 10/5/2022, eq (2)
matrix_tt = matrix(1:2,1:2);
matrix_t = matrix(1:2,3);
matrix_z = matrix(3,1:2).'; % HH transposes this so _t and _z are both column vectors
matrix_zz = matrix(3,3);
end

