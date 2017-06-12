function [agOutputVertexList] = takeNeighbours(agVertexList, agIndexList, cmpHds, nNeighbours)
%TAKENEIGHBOURS Takes only the nNeighbours from the list of agVertexList
%and outputs them, based on the cmpHds values.
%Input:
%  cmpHds: Compared heat descriptors of vertices, and the norm of that
%  value

%% init values
nElems = size(agVertexList,1);
dataPoints = zeros(nElems,1);
for i = 1:nElems
    idx = agIndexList(i);
    dataPoints(i) = cmpHds(idx);
end

%% find the minimum nNeighbours
maxVal = max(dataPoints) * 20;
minIndexPositions = [];

for i = 1:nNeighbours
    [~,minPos] = min(dataPoints);
    minIndexPositions = [minIndexPositions; minPos];
    dataPoints(minPos) = maxVal;
end

%% set the outputVertexList with the required ones
agOutputVertexList = agVertexList(minIndexPositions,:);

end