function [vertices, avgDist] = addSuccessiveNoise(vertices, percentageNoise)
%ADDNOISETOTRAJECTORYPOSES Adds (x,y) noise to vertices, it is done by
%percentage amount of the average trajectory step. The noise is accumulated
%over each step.
% Input:
%  vertices: array of struct for each pose information
%  percentageNoise: 0.05 for 5% pertubation in the pose information


vCount = size(vertices,2);
totDistn = 0;
poseDiff = zeros(vCount-1,2);
% -- calculate the average pose distance
for i = 2:vCount
    x1 = vertices(i-1).x;
    x2 = vertices(i).x;
    y1 = vertices(i-1).y;
    y2 = vertices(i).y;
    % -- calculate the pose distance
    poseDiff(i-1,1) = x2 - x1;
    poseDiff(i-1,2) = y2 - y1;
end
poseDiffsq = poseDiff.^2;
poseDiffsqSum = sum(poseDiffsq')';
poseDiffsqSum = poseDiffsqSum.^0.5;
avgDist = mean(poseDiffsqSum);

noiseLevel = percentageNoise * avgDist;
noiseData_x = (rand(1,vCount) - 0.5)*2*noiseLevel;
noiseData_y = (rand(1,vCount) - 0.5)*2*noiseLevel;
% Note: rand generates random value between 0 and 1, so subtracting 0.5
% would generate the random value between -0.5 to +0.5. And then
% multiplying it with 2 and noiseLevel would give the data we need.

for i = 2:vCount
    vertices(i).x = vertices(i-1).x + noiseData_x(i) + poseDiff(i-1,1);
    vertices(i).y = vertices(i-1).y + noiseData_y(i) + poseDiff(i-1,2);
end
end