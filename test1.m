% check if the plotting function works well.

% Detailed:
% In this test, we load up the two mat files corresponding to the two
% graphs, and then find out the ground truth associations. After they have
% been found, we then consider to plot each of the 10th association in the
% figure. This shows if the plotting function works well, and there are no
% glitches in the plotting functions.

clc;

%% get sparse landmark association

file1 = 'input/databaseFile-1.mat';
file2 = 'input/databaseFile-2.mat';
landmarkAssc = findLandmarkAssociations(file1, file2);
fprintf(1,'Status: Found landmark Associations\n');
poseAssc = findPoseAssociations(file1, file2);
fprintf(1,'Status: Found pose Associations\n');


%% file 1
load(file1);
vertices1 = vertices;
landmarks1 = landmarks;
edges1 = edges;
ledges1 = ledges;
fprintf(1,'Status: Loaded file1 data\n');

%% file 2

load(file2);
vertices2 = vertices;
landmarks2 = landmarks;
edges2 = edges;
ledges2 = ledges;
fprintf(1,'Status: Loaded file2 data\n');

%% visualize data

% sparsify the denseAssc;
denseAssc = poseAssc;
for i = 1:length(denseAssc)
    if (mod(i,10) ~= 0)
        denseAssc(i) = 0;
    end
end


figure('name','data-Visualization');
lA1 = landmarkAssc > 0;
lA2tmp = landmarkAssc(landmarkAssc>0);
lA2 = zeros(size(landmarks2,2));
lA2(lA2tmp) = 1;
depth_z = 50;
x_offset = 500;
y_offset = 0;
plotGraph(vertices1, landmarks1, lA1, edges1, ledges1, 0, 0, depth_z, 'blue','red');
plotGraph(vertices2, landmarks2, lA2, edges2, ledges2, x_offset, y_offset, depth_z, 'green','cyan');
plotDenseMap(denseAssc, vertices1, vertices2, landmarks1, landmarks2, x_offset,y_offset, depth_z);
