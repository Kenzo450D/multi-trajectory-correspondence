function [vertexNeighbours] = getNeighbours(leData, lData,landmarkIdx)

% -- get the landmark id
lId = lData(1,landmarkIdx);

% -- get the vertexPoints
row1Points = leData(2,ismember(leData(1,:),lId));
row2Points = leData(1,ismember(leData(2,:),lId));

% -- add both results
vertexNeighbours = [row1Points, row2Points];
vertexNeighbours = unique(vertexNeighbours);

end