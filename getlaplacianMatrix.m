function [L] = getlaplacianMatrix(A)
%GETLAPLACIANMATRIX Returns the laplacian matrix from Adjacency Matrix

D = sum(A);
L = diag(D) - A;

end
