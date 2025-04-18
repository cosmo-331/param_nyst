"main.m" is the main script that we want to run. 

"param_nyst.m" is the function that computes a list of indices for all parameter values, like AdaCUR. 

"RPCholesky.m" is my implementation of RPCholesky. 

"RPCholesky_preindex.m" is the RPCholesky algorithm when we already have an index set and want to append new indices to it. 

"reduce_column.m", "reduce_column_chol.m", and "reduce_column_sketch.m" are three different implementations in which we remove redundant indices. "reduce_column_sketch" is recommended.

"ACA.m" and "param_ACA.m" are my implementations of the algorithms described in Kressner's paper https://arxiv.org/abs/2001.09187.

"insertcolumn.m" is downloaded from MATLAB file exchange https://uk.mathworks.com/matlabcentral/fileexchange/90835-updating-thin-qr-factorization.
