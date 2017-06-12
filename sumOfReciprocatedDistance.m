function [vlDist] = sumOfReciprocatedDistance(vData,lData, eData, leData)

% here distance is the average size of odometry step.
% we consider that a landmark won't have more than 100 vertices as
% neighbours
%% init

maxNeighbours = 100;

%% calculate the average step size in eData

avgStepSize = getAverageStepSize(eData);

%% get the neighbours
vCount = size(vData,2);
lCount = size(lData,2);

landmarkIdx = (1:length(lData));
neighbourList  = zeros(lCount,maxNeighbours);
for i = 1:lCount
    landmarkIdx = i;
    nbrs = getNeighbours(leData, lData,landmarkIdx);
    neighbourList(i,1:length(nbrs)) = nbrs;
end

%% form matrix to get steps from vertex to lData

stepsVertexLandmark = zeros(vCount,lCount);
for i = 1:vCount
    % find the landmarks which are neighbours
    row1Neighbours = leData(2,ismember(leData(1,:),vData(1,i)));
    row2Neighbours = leData(1,ismember(leData(2,:),vData(1,i)));
    landmarkNeighbours = [row1Neighbours, row2Neighbours];
    for j = 1:lCount
        if (ismember(landmarkNeighbours,j))
            stepsVertexLandmark(i,j) = 1;
        else
            stepsVertexLandmark(i,j) = abs(vData(1,i) - lData(1,j));
        end
    end
end

