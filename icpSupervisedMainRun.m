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

%% add landmark noise
tmp              = landmarkAssc(2);
landmarkAssc(2)  = landmarkAssc(40);
landmarkAssc(40) = tmp;
tmp              = landmarkAssc(10);
landmarkAssc(10) = landmarkAssc(20);
landmarkAssc(20) = tmp;
tmp              = landmarkAssc(15);
landmarkAssc(15) = landmarkAssc(25);
landmarkAssc(25) = tmp;

%% add pose noise
[vertices1] = addSuccessiveNoise(vertices1,1);
[vertices2] = addSuccessiveNoise(vertices2,1);

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

%% run ICP on pose space

[Ricp, Ticp, ER, t, match] = icp(poseData1,poseData2,15);

%%  convert match to form of matchesMade
tj2v1 = 1;
tj2v2 = size(vData2,2);
tj1v1 = 1;
tj1v2 = size(vData1,2);
newMatchesMade = zeros(size(match,2),2);
newMatchesMade(:,2) = 1:(tj2v2 - tj2v1 + 1);
newMatchesMade(:,2) = newMatchesMade(:,2) + tj2v1 - 1; % the '1' is to adjust for the matlab indexing (starts from 1 and not 0)
newMatchesMade(:,1) = match + tj1v1 - 1;

%% visualize data
%{
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
%}
%% plot

totalVertexCount1 = size(vertices1,2) + size(landmarks1,2);
denseAssc = zeros(totalVertexCount1,1);
if (~isempty(newMatchesMade))
    denseAssc(vData1(1,newMatchesMade(:,1))) = vData2(1,newMatchesMade(:,2));
end
% 
% % -- get the poses which has an association
% poseIdx = gtposeAssc>0;
% posesToCompare = gtposeAssc(poseIdx);

% to match the posesToCompare to that of the newMatchesMade, we could crop
% the output by newMatchesMade to the region that we require, and take the
% result of only that section

% -- map the posesToCompare to poseGraph2
% matchesSelected = newMatchesMade(1);

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

% -- calculate ATE and precision
[ate] = calculateATE(vData1, vData2, newMatchesMade);
[p, relaxPrecision] = calculatePrecision(gtMatchesMade, newMatchesMade);

outFile = 'atePrecision.txt';
fid = fopen(outFile,'a')';

fprintf(1,'%d %f %f\n',ate,p,relaxPrecision);
fprintf(fid,'%f %f %f\n',ate,p,relaxPrecision);
fclose(fid);