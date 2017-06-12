function [distn, vertexDescriptors, t_scale] = heatEmbeddingDescriptor(L, t_scale_indicator, t_scale, n_perc_eigval_indx, vertices, landmarks, edges, ledges, landmarkAssc)
%HEATEMBEDDINGDESCRIPTOR Heat dissipation distance for several t_scale params
% L        : Laplacian Matrix
% t_scale  : vector of time scale parameters for heat kernel embedding
% vertices : vertices of the pose grap
% edges    : edges of the pose graph
% landmarks: landmarks of the pose graph
% ledges   : landmark edges of the pose graph
% -------------------------------------------------------------------------
% distn    : A matrix of ntScales,totalVertexCount,totalVertexCount
% vertexDescriptors: A matrix of ntScales x commonLandmarks for every
%                    vertex

%% initialization
vCount    = size(vertices,2);
eCount    = size(edges,2);
lCount    = size(landmarks,2);
nlandAssc = nnz(landmarkAssc);

%% eigen decomposition of the the matrices
[eigenVectors, eigenValues] = eig(L);
eigenValues                 = diag(eigenValues);
nEigenValues                = length(eigenValues);

%% as well some heat embedding
if (isempty(t_scale))
    t_scale=-log(t_scale_indicator)/eigenValues(ceil(size(L,1)*n_perc_eigval_indx));
end
ntScales = length(t_scale);
% -- plot the scales
% figure('name','scale diffusion');
% for i = 1:ntScales
% plot(eigenValues,exp(-eigenValues*t_scale(i)));
% hold on;
% end;
% hold off;
% -- end plot
totalVertexCount = vCount + lCount;
distn = zeros(ntScales,totalVertexCount,totalVertexCount);
for i=1:ntScales
    [ heatDistance ] = distanceEmbedding( eigenVectors,eigenValues,nEigenValues, t_scale(i));
    distn(i,:,:) = heatDistance;
end

% -- we want a matrix for every vertex
% the matrix would be sized ntScales x commonLandmarks

vertexDescriptors = zeros(ntScales, nlandAssc, vCount);
for i = 1:vCount
    % -- set each of vertices
    for j = 1:ntScales
        for k = 1:nlandAssc
            vertexDescriptors(j,k,i) = distn(j,vertices(i).id,landmarks(landmarkAssc(k)).id);
        end
    end
end

end