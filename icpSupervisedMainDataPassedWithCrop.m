function [newMatchesMade] = icpSupervisedMainDataPassedWithCrop(vertices1, vertices2, landmarks1, landmarks2, tj1v1, tj1v2, tj2v1, tj2v2, landmarkAssc, atePrecisionFile, outMatICPFileName)
%ICPMAINDATAPASSED It is a data passed version of ICP Only code. This
%accepts the vertex and edge information, and runs unsupervised ICP on it. 

% -- convert to matrix
vData1 = getVertexMatrixForm(vertices1);
vData2 = getVertexMatrixForm(vertices2);

%% initialize the icp with the landmarkAssc
[R, T] = icpForLandmarkAssc(landmarkAssc, landmarks1, landmarks2);

%% transform the point cloud for vData2

% -- initialize poseData1 and poseData2
poseData1 = vData1(2:3,:);
poseData2 = vData2(2:3,:);
% -- zero Padding
zp = zeros(1,size(vData1, 2));
poseData1 = [poseData1;zp];
zp = zeros(1,size(vData2, 2));
poseData2 = [poseData2;zp];

% -- transform poseData2
poseData2 = R * poseData2 + repmat(T, 1, size(poseData2,2));

%% truncate the poseData to the range declared by start and end points
% [maxVal, maxIdx] = max(gtPoseAssc);
% tj2v1 = gtPoseAssc(1);
% tj2v2 = maxVal;
% tj1v1 = 1;
% tj1v2 = maxIdx;
poseData1 = poseData1(:,tj1v1:tj1v2);
poseData2 = poseData2(:,tj2v1:tj2v2);


%% run ICP on pose space

[Ricp, Ticp, ER, t, match] = icp(poseData1,poseData2,15);

%%  convert match to form of matchesMade

newMatchesMade = zeros(size(match,2),2);
newMatchesMade(:,2) = 1:(tj2v2 - tj2v1 + 1);
newMatchesMade(:,2) = newMatchesMade(:,2) + tj2v1 - 1; % the '1' is to adjust for the matlab indexing (starts from 1 and not 0)
newMatchesMade(:,1) = match + tj1v1 - 1;

%% plot
%{
totalVertexCount1 = size(vertices1,2) + size(landmarks1,2);
denseAssc = zeros(totalVertexCount1,1);
if (~isempty(newMatchesMade))
    denseAssc(vData1(1,newMatchesMade(:,1))) = vData2(1,newMatchesMade(:,2));
end

% -- plot the data
figure('name','match plot');
lA1 = landmarkAssc > 0;
lA2tmp = landmarkAssc(landmarkAssc>0);
lA2 = zeros(size(landmarks2,2));
lA2(lA2tmp) = 1;
depth_z = 50;
x_offset = 500;
y_offset = 0;
color1 = jet(size(vertices1,2));
plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, color1,'red');
color2 = getColor2(color1, size(vertices2,2), denseAssc);
% plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, 'blue','red');
plotGraph(vertices2, landmarks2, lA2, x_offset, y_offset, depth_z, color2,'cyan');
plotDenseMap(denseAssc, vertices1, vertices2, landmarks1, landmarks2, x_offset,y_offset, depth_z,5);
%}

%% ate and precision

% -- convert gtPoseAssc to gtMatchesMade
% poseIdx = gtPoseAssc>0;
% posesToCompare = gtPoseAssc(poseIdx);
% gtMatchesMade = zeros(length(posesToCompare),2);
% for i = 1:length(posesToCompare)
%     if (poseIdx(i) ~= 0)
%         gtMatchesMade(i,1) = i;
%         gtMatchesMade(i,2) = gtPoseAssc(i);
%     end
% end
% 
% 
% [ate] = calculateATE(vData1, vData2, newMatchesMade);
% [p, relaxPrecision] = calculatePrecision(gtMatchesMade, newMatchesMade, (size(vData1,2)*0.001));
% 
% outFile = atePrecisionFile;
% fid = fopen(outFile,'a')';
% 
% fprintf(1,'ICP: %f %f %f\n',ate,p,relaxPrecision);
% fprintf(fid,'%f %f %f\n',ate,p,relaxPrecision);
% fclose(fid);
% 
% %% save the matches made in the file
% icpMatchesMade = newMatchesMade;
% save(outMatICPFileName,'icpMatchesMade');

end