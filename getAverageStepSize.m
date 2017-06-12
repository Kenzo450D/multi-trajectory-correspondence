function [avgStepSize] = getAverageStepSize(eData)
%GETAVERAGESTEPSIZE Returns the average step size in a trajectory
% Input:
%   eData: matrix of all the edges
% Output:
%   avgStepSize: average Step Size of the trajectory
eCount = size(eData,2);
totalStepSize = 0;
totalStepCount = 0;
for i = 1:eCount
    %check for odometry edge (usually if consecutive or with one gap)
    v1 = eData(1,i);
    v2 = eData(2,i);
    poseIdxDiff = v2 - v1;
    if (poseIdxDiff == 1 || poseIdxDiff == 0)
        % is an odometry edge
        stepsize = sqrt(eData(3,i)*eData(3,i) + eData(4,i) * eData(4,i));
        totalStepSize = totalStepSize + stepsize;
        totalStepCount = totalStepCount + 1;
    end
end

avgStepSize = totalStepSize / totalStepCount;