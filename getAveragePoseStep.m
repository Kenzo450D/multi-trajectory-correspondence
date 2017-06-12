function [avgPoseStep] = getAveragePoseStep(vData)
%GETAVERAGEPOSESTEP Returns the average pose step of all odometry edges

vCount = size(vData,2);
avgPoseStep = 0;
for i = 2:vCount
    px2 = vData(2,i);
    py2 = vData(3,i);
    px1 = vData(2,(i-1));
    py1 = vData(3,(i-1));
    diffX = px2 - px1;
    diffY = py2 - py1;
    stepSize = sqrt(diffX*diffX + diffY*diffY);
    avgPoseStep = avgPoseStep + stepSize;
end
avgPoseStep = avgPoseStep / (vCount - 1);
end