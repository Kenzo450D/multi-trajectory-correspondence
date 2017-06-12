clear; clc;


%% -- Victoria park
file1 = 'input/databaseFile-1.mat';
file2 = 'input/databaseFile-2.mat';
optimisedFile1 = 'input/optimised/victoriaParkOpDataset-1.mat';
optimisedFile2 = 'input/optimised/victoriaParkOpDataset-2.mat';


%% -- load unoptimised file
% ---- load file 1
load(file1);
vertices1 = vertices;
landmarks1 = landmarks;
edges1 = edges;
ledges1 = ledges;
% ---- load file 2
load(file2);
vertices2 = vertices;
landmarks2 = landmarks;
edges2 = edges;
ledges2 = ledges;



%% find gt associations
landmarkAssc = findLandmarkAssociations(file1, file2);
[gtposeAssc, gtDenseAssc] = findPoseAssociations(file1, file2);

%% visualize data
% figure('name','compare plot');
% lA1 = landmarkAssc > 0;
% lA2tmp = landmarkAssc(landmarkAssc>0);
% lA2 = zeros(size(landmarks2,2));
% lA2(lA2tmp) = 1;
% depth_z = 50;
% x_offset = 500;
% y_offset = 0;
% color1 = jet(size(vertices1,2));
% plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, color1,'red');
% color2 = getColor2(color1, size(vertices2,2), gtDenseAssc);
% % plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, 'blue','red');
% plotGraph(vertices2, landmarks2, lA2, x_offset, y_offset, depth_z, color2,'cyan');
% plotDenseMap(gtDenseAssc, vertices1, vertices2, landmarks1, landmarks2, x_offset,y_offset, depth_z);
% return;

%% call the ICP function
%[vertices1] = addNoiseToTrajectoryPoses(vertices1,10);
%[vertices2] = addNoiseToTrajectoryPoses(vertices2,10);

[Ricp, Ticp, ER, t, newMatchesMade] = runICPOnly(vertices1, vertices2);

% -- convert to matrix
[vData1, lData1, edges1, ledges1] = getMatrixForm(vertices1, landmarks1, edges1, ledges1);
[vData2, lData2, edges2, ledges2] = getMatrixForm(vertices2, landmarks2, edges2, ledges2);


totalVertexCount1 = size(vertices1,2) + size(landmarks1,2);
denseAssc = zeros(totalVertexCount1,1);
if (~isempty(newMatchesMade))
    denseAssc(vData1(1,newMatchesMade(:,1))) = vData2(1,newMatchesMade(:,2));
end

% -- get the poses which has an association
poseIdx = gtposeAssc>0;
posesToCompare = gtposeAssc(poseIdx);

% to match the posesToCompare to that of the newMatchesMade, we could crop
% the output by newMatchesMade to the region that we require, and take the
% result of only that section

% -- map the posesToCompare to poseGraph2
matchesSelected = newMatchesMade(1);

% -- plot the data
% figure('name','match plot');
% lA1 = landmarkAssc > 0;
% lA2tmp = landmarkAssc(landmarkAssc>0);
% lA2 = zeros(size(landmarks2,2));
% lA2(lA2tmp) = 1;
% depth_z = 50;
% x_offset = 500;
% y_offset = 0;
% color1 = jet(size(vertices1,2));
% plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, color1,'red');
% color2 = getColor2(color1, size(vertices2,2), denseAssc);
% % plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, 'blue','red');
% plotGraph(vertices2, landmarks2, lA2, x_offset, y_offset, depth_z, color2,'cyan');
% plotDenseMap(denseAssc, vertices1, vertices2, landmarks1, landmarks2, x_offset,y_offset, depth_z);
% return;

%% ate and precision

% -- convert gtPoseAssc to gtMatchesMade
poseIdx = gtposeAssc>0;
posesToCompare = gtposeAssc(poseIdx);
gtMatchesMade = zeros(length(posesToCompare),2);
for i = 1:length(posesToCompare)
    if (poseIdx(i) ~= 0)
        gtMatchesMade(i,1) = i;
        gtMatchesMade(i,2) = gtposeAssc(i);
    end
end

[ate] = calculateATE(vData1, vData2, newMatchesMade);
[p, relaxPrecision] = calculatePrecision(gtMatchesMade, newMatchesMade);

outFile = 'atePrecision.txt';
fid = fopen(outFile,'a')';

fprintf(1,'%d %d %d\n',ate,p,relaxPrecision);
fprintf(fid,'%d %d %d\n',ate,p,relaxPrecision);
fclose(fid);