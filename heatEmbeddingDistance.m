function [distn] = heatEmbeddingDistance(A, t_scale, vertices, landmarks, edges, ledges)
%HEATEMBEDDINGDISTANCE Heat dissipation distance for several t_scale params
% A: incidence matrix
% t_scale: vector of time scale parameters for heat kernel embedding
% edges: edges of the pose graph
% vertices: vertices of the pose graph

% -- construct laplacian matrix
Asp = sparse(A);
L = Asp*Asp';
L = full(L);

% -- eigen decomposition of the the matrices
[eigenVectors, eigenValues]=eig(L);
eigenValues=diag(eigenValues);
nEigenValues=length(eigenValues);

% -- as well some heat embedding
ntScales = length(t_scale);
totalVertexCount = size(vertices,2) + size(landmarks,2);
distn = zeros(ntScales,totalVertexCount,totalVertexCount);
for i=1:ntScales
    [ heatDistance ] = distanceEmbedding( eigenVectors,eigenValues,nEigenValues, edges, ledges, t_scale(i), vertices, landmarks);
    distn(i,:,:) = heatDistance;
end
end
