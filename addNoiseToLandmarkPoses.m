function [landmarks] = addNoiseToLandmarkPoses(landmarks, percentageNoise, avgDist)
%ADDNOISETOTRAJECTORYPOSES Adds (x,y) noise to vertices, it is done by
%percentage amount of the average trajectory step
% Input:
%  vertices: array of struct for each pose information
%  percentageNoise: 0.05 for 5% pertubation in the pose information
%  avgDist: average Odometry distance


lCount = size(landmarks,2);
totDistn = 0;
% -- calculate the average pose distance
for i = 2:lCount
    x1 = landmarks(i-1).x;
    x2 = landmarks(i).x;
    y1 = landmarks(i-1).y;
    y2 = landmarks(i).y;
    % -- calculate the pose distance
    distn = sqrt((x1 - x2)^2 + (y1 - y2)^2);
    totDistn = totDistn + distn;
end
avgDist = totDistn / (lCount - 1);

noiseLevel = percentageNoise * avgDist;
noiseData_x = (rand(1,lCount) - 0.5)*2*noiseLevel;
noiseData_y = (rand(1,lCount) - 0.5)*2*noiseLevel;
% Note: rand generates random value between 0 and 1, so subtracting 0.5
% would generate the random value between -0.5 to +0.5. And then
% multiplying it with 2 and noiseLevel would give the data we need.

for i = 1:lCount
    landmarks(i).x = landmarks(i).x + noiseData_x(i);
    landmarks(i).y = landmarks(i).y + noiseData_y(i);
end
end