function makeDataset( vertices, landmarks, edges, ledges, outputFileBase )
%MAKEDATASET Deletes the first and last 1500 vertices, and makes two mat
%files out of them

%% make the outputFile names
outputFile1 = [outputFileBase,'-1.mat'];
outputFile2 = [outputFileBase,'-2.mat'];
outputG2oFile1 = [outputFileBase,'-1.g2o'];
outputG2oFile2 = [outputFileBase,'-2.g2o'];

%% convert the data to matrices
% -- get the size of each of the inputs
vCount = size(vertices,2);
lCount = size(landmarks,2);
eCount = size(edges,2);
leCount = size(ledges,2);


% -- convert the structure data to an array
[vertexData, landmarkData, edgeData, landmarkEdgeData] = getMatrixForm(vertices, landmarks, edges, ledges);

%% take the backup of the data
initVertexData = vertexData;
initEdgeData = edgeData;
initLandmarkData = landmarkData;
initLandmarkEdgeData = landmarkEdgeData;

%% delete the first 1500 vertices
% -- vertices <= 1500
% ---- find the indexes
vIdx = find(vertexData(1,:) <= 1500);
% ---- remove the vertices
vertexData(:,vIdx) = [];

% -- landmarks <= 1500
% ---- find the indexes
vIdx = find(landmarkData(1,:) <= 1500);
% ---- remove the landmarks
landmarkData(:,vIdx) = [];

% -- edges with vertices <= 1500
% ---- find the edges with first vertex <= 1500
eIdx = find(edgeData(1,:) <= 1500);
% ---- remove the edges
edgeData(:,eIdx) = [];
% ---- find the edges with second vertex <= 1500
eIdx = find(edgeData(2,:) <= 1500);
% ---- remove the edges
edgeData(:,eIdx) = [];

% -- landmark edges with the vertices <= 1500
% ---- find the landmark edges with first vertex <= 1500
eIdx = find(landmarkEdgeData(1,:) <= 1500);
% ---- remove the landmark edges
landmarkEdgeData(:,eIdx) = [];
% ---- find the landmark edges with second vertex <= 1500
eIdx = find(landmarkEdgeData(2,:) <= 1500);
% ---- remove the landmark edges
landmarkEdgeData(:,eIdx) = [];

% -- get the new set of indexes
[vertexData, landmarkData, edgeData, landmarkEdgeData] = fixIndexes(vertexData, landmarkData, edgeData, landmarkEdgeData);

fprintf(1,'vCount: %d\n', size(vertexData,2));

% -- convert the data to mat file, and save it
convertToMatFile(vertexData, landmarkData, edgeData, landmarkEdgeData, outputFile1, outputG2oFile1);

%% delete the last 1500 vertices
vertexData = initVertexData;
edgeData = initEdgeData;
landmarkData = initLandmarkData;
landmarkEdgeData = initLandmarkEdgeData;


% -- find the index
maxIdx = max(vertexData(1,:));
if ( maxIdx > max(landmarkData(1,:)))
    maxIdx = max(landmarkData(1,:));
end
limitIdx = maxIdx - 1500;

% -- delete the not required vertices
vIdx = find(vertexData(1,:) > limitIdx);
vertexData(:,vIdx) = [];
vIdx = find(landmarkData(1,:) > limitIdx);
landmarkData(:,vIdx) = [];
eIdx = find(edgeData(1,:) > limitIdx);
edgeData(:,eIdx) = [];
eIdx = find(edgeData(2,:) > limitIdx);
edgeData(:,eIdx) = [];
eIdx = find(landmarkEdgeData(1,:) > limitIdx);
landmarkEdgeData(:,eIdx) = [];
eIdx = find(landmarkEdgeData(2,:) > limitIdx);
landmarkEdgeData(:,eIdx) = [];

% -- get the new set of indexes
[vertexData, landmarkData, edgeData, landmarkEdgeData] = fixIndexes(vertexData, landmarkData, edgeData, landmarkEdgeData);

% -- convert the data to mat file, and save it
convertToMatFile(vertexData, landmarkData, edgeData, landmarkEdgeData, outputFile2, outputG2oFile2);



end

