% check if the vertex Descriptors are working fine

% Detailed:
% In this test, we load up the two mat files corresponding to the two
% graphs, and then find out the ground truth associations. After they have
% been found, we then find out the vertex descriptors of all of the
% vertices of the two graphs. Once the vertex descriptors have been found
% out, we compare between them, and find out whether they match up for the
% ground truth neighbours. We plot 5 cases:
% 1:2 matches which are close to landmarks
% 3 non-match while are close to landmarks
% 4 match which are not close to landmark
% 5 non-match which aren't close to landmark

clc;

%% get sparse landmark association

file1 = 'input/databaseFile-1.mat';
file2 = 'input/databaseFile-2.mat';
landmarkAssc = findLandmarkAssociations(file1, file2);
fprintf(1,'Status: Found landmark Associations\n');
poseAssc = findPoseAssociations(file1, file2);
fprintf(1,'Status: Found pose Associations\n');

%% initialize with g2o optimised maps

%t_scale = [0.001,0.01,1,10,50];
t_scale_indicator = [0.99,0.95,0.92,0.9,0.85, 0.8];
n_perc_eigval_indx = 0.5/100;
%t_scale = [0.001,0.01,1,10,100];

%% file 1
load(file1);
% ---- make the adjacency matrix
A = makeAdjacencyMatrix(vertices, landmarks, edges, ledges);
% ---- get the laplacian matrix
L = getlaplacianMatrix(A);

% ---- prepare associations
% find the poses in map1 which have an association with the other map
l1Assc = landmarkAssc>0;
tmpVar = [1:size(landmarks,2)];
l1Assc = tmpVar(l1Assc);

% ---- prepare heat embedding descriptor
%vd1 : vertexDescriptors1
[distn1, vd1] = heatEmbeddingDescriptor(L, t_scale_indicator, n_perc_eigval_indx, vertices, landmarks, edges, ledges, l1Assc);
vertices1 = vertices;
landmarks1 = landmarks;
edges1 = edges;
ledges1 = ledges;
fprintf(1,'Status: Loaded file1 data\n');


%% file 2

load(file2);
A = makeAdjacencyMatrix(vertices, landmarks, edges, ledges);
L = getlaplacianMatrix(A);

% -- prepare associations
% find the poses in map2 which have an association in map1
l2Assc = landmarkAssc>0;
l2Assc = landmarkAssc(l2Assc); % TODO: Check this

[distn2,vd2] = heatEmbeddingDescriptor(L, t_scale_indicator, n_perc_eigval_indx, vertices, landmarks, edges, ledges, l2Assc);
vertices2 = vertices;
landmarks2 = landmarks;
edges2 = edges;
ledges2 = ledges;
fprintf(1,'Status: Loaded file2 data\n');

%% get data into matrices
vCount1 = size(vertices1, 2);
vCount2 = size(vertices2, 2);
lCount1 = size(landmarks1, 2);
lCount2 = size(landmarks2, 2);
totalVertexCount1 = vCount1 + lCount1;
totalVertexCount2 = vCount2 + lCount2;
[vData1, lData1, eData1, leData1] = getMatrixForm(vertices1, landmarks1, edges1, ledges1);
[vData2, lData2, eData2, leData2] = getMatrixForm(vertices2, landmarks2, edges2, ledges2);

%% create reverse vertex map
vMap1 = reverseVertexMap(vData1);
vMap2 = reverseVertexMap(vData2);


%% get neighbours of both files

landmarkIdx = (1:length(lData1));
[vn1] = getNeighbours(leData1, lData1,landmarkIdx);
landmarkIdx = (1:length(lData2));
[vn2] = getNeighbours(leData2, lData2, landmarkIdx);
vrtxCombo = zeros(5,2);
vn1Count = size(vn1);
vn2Count = size(vn2);

% consider two neighbours close to landmarks
vrtxCombo(1,1) = 3110;
vrtxCombo(2,1) = 3115;
vrtxCombo(1,2) = vData2(1,poseAssc(vMap1(vrtxCombo(1,1))));
vrtxCombo(2,2) = vData2(1,poseAssc(vMap1(vrtxCombo(2,1))));

% consider two neighbours not close to landmarks
vrtxCombo(3,1) = 600;
vrtxCombo(4,1) = 700;
vrtxCombo(3,2) = vData2(1,poseAssc(vMap1(vrtxCombo(3,1))));
vrtxCombo(4,2) = vData2(1,poseAssc(vMap1(vrtxCombo(4,1))));

% consider a pair not related to one another (close to landmarks)
vrtxCombo(5,1) = 3112;
vrtxCombo(5,2) = vData2(1,poseAssc(vMap1(vrtxCombo(5,1)))) + 3;


%% plot the vertexDescriptors of them
for i = 1:5
    str = sprintf('%d: Vertex Descriptor of %d in Graph1, %d in Graph2',i,vrtxCombo(i,1), vrtxCombo(i,2));
    figure('name',str);
    showVd1 = vd1(:,:,vMap1(vrtxCombo(i,1)));
    showVd2 = vd2(:,:,vMap2(vrtxCombo(i,2)));
    subplot(1,2,1);
    imagesc(showVd1);
    subplot(1,2,2);
    imagesc(showVd2);
end