function [vertices, avgDist] = addNoiseToTrajectoryPoses(vertices, percentageNoise)
%ADDNOISETOTRAJECTORYPOSES Adds (x,y) noise to vertices, it is done by
%percentage amount of the average trajectory step
% Input:
%  vertices: array of struct for each pose information
%  percentageNoise: 0.05 for 5% pertubation in the pose information


vCount = size(vertices,2);
totDistn = 0;
% -- calculate the average pose distance
for i = 2:vCount
    x1 = vertices(i-1).x;
    x2 = vertices(i).x;
    y1 = vertices(i-1).y;
    y2 = vertices(i).y;
    % -- calculate the pose distance
    distn = sqrt((x1 - x2)^2 + (y1 - y2)^2);
    totDistn = totDistn + distn;
end
avgDist = totDistn / (vCount - 1);

noiseLevel = percentageNoise * avgDist;
noiseData_x = (rand(1,vCount) - 0.5)*2*noiseLevel;
noiseData_y = (rand(1,vCount) - 0.5)*2*noiseLevel;
% Note: rand generates random value between 0 and 1, so subtracting 0.5
% would generate the random value between -0.5 to +0.5. And then
% multiplying it with 2 and noiseLevel would give the data we need.

for i = 1:vCount
    vertices(i).x = vertices(i).x + noiseData_x(i);
    vertices(i).y = vertices(i).y + noiseData_y(i);
end
end