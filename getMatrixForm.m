function [vertexData, landmarkData, edgeData, landmarkEdgeData] = getMatrixForm(vertices, landmarks, edges, ledges)

vCount = size(vertices,2);
lCount = size(landmarks,2);
eCount = size(edges,2);
leCount = size(ledges,2);

vertexData = zeros(4,vCount);
for i = 1:vCount
    vertexData(1,i) = vertices(i).id;
    vertexData(2,i) = vertices(i).x;
    vertexData(3,i) = vertices(i).y;
    vertexData(4,i) = vertices(i).o;
end
landmarkData = zeros(3,lCount);
for i = 1:lCount
    landmarkData(1,i) = landmarks(i).id;
    landmarkData(2,i) = landmarks(i).x;
    landmarkData(3,i) = landmarks(i).y;
end
edgeData = zeros(11,eCount);
for i = 1:eCount
    edgeData(1,i) = edges(i).v1;
    edgeData(2,i) = edges(i).v2;
    edgeData(3,i) = edges(i).dx;
    edgeData(4,i) = edges(i).dy;
    edgeData(5,i) = edges(i).dth;
    %fill Covariance Matrix
    edgeData(6,i)  = edges(i).covMatrix(1,1);
    edgeData(7,i)  = edges(i).covMatrix(1,2);
    edgeData(8,i)  = edges(i).covMatrix(1,3);
    edgeData(9,i)  = edges(i).covMatrix(2,2);
    edgeData(10,i) = edges(i).covMatrix(3,2);
    edgeData(11,i) = edges(i).covMatrix(3,3);
end
landmarkEdgeData = zeros(7,leCount);
for i = 1:1:leCount
    landmarkEdgeData(1,i) = ledges(i).v1;
    landmarkEdgeData(2,i) = ledges(i).v2;
    landmarkEdgeData(3,i) = ledges(i).dx;
    landmarkEdgeData(4,i) = ledges(i).dy;
    landmarkEdgeData(5,i) = ledges(i).covMatrix(1,1);
    landmarkEdgeData(6,i) = ledges(i).covMatrix(1,2);
    landmarkEdgeData(7,i) = ledges(i).covMatrix(2,2);
end


end