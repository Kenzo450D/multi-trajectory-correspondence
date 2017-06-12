function [vertexData, landmarkData, edgeData, landmarkEdgeData] = fixIndexes(vertexData, landmarkData, edgeData, landmarkEdgeData)

% -- find the minimum of the vertices, truncate the values accordingly
minVertexId = min(vertexData(1,:));
if (minVertexId > min(landmarkData(1,:)))
    minVertexId = min(landmarkData(1,:));
end
% -- matlab vertices start from 1
minVertexId = minVertexId - 1; 

% -- adjust the rest of the edges and vertices accordingly so that the
% minimum is zero

vertexData(1,:) = vertexData(1,:) - minVertexId;
landmarkData(1,:) = landmarkData(1,:) - minVertexId;
edgeData(1,:) = edgeData(1,:) - minVertexId;
edgeData(2,:) = edgeData(2,:) - minVertexId;
landmarkEdgeData(1,:) = landmarkEdgeData(1,:) - minVertexId;
landmarkEdgeData(2,:) = landmarkEdgeData(2,:) - minVertexId;
end