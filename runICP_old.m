clear;
clc;

%% get sparse landmark association

file1 = 'input/databaseFile-1.mat';
file2 = 'input/databaseFile-2.mat';
landmarkAssc = findLandmarkAssociations(file1, file2);
[poseAssc, gtDenseAssc] = findPoseAssociations(file1, file2);

tic;
%% load each file to find 
% =========================================================================
% FILE 1
% =========================================================================
% -- load the file
load(file1);
% -- add noise to each of the files
[vertices] = addNoiseToTrajectoryPoses(vertices,10);

vertices1 = vertices;
landmarks1 = landmarks;
edges1 = edges;
ledges1 = ledges;

% =========================================================================
% FILE 2
% =========================================================================
load(file2);
% -- add noise to each of the files
[vertices] = addNoiseToTrajectoryPoses(vertices,10);

vertices2 = vertices;
landmarks2 = landmarks;
edges2 = edges;
ledges2 = ledges;
toc;

%% convert vertices to matrix data
[vData1] = getVertexMatrixForm(vertices1);
[vData2] = getVertexMatrixForm(vertices2);

%% convert matrix data to 3 dim
vData1 = vData1(2:3,:);
vData2 = vData2(2:3,:);
o = zeros(1,size(vData1,2)); % last row
vData1 = [vData1; o];
o = zeros(1,size(vData2,2)); % last row
vData2 = [vData2; o];

[Ricp Ticp ER t] = icp(vData1, vData2, 150);


%% quantitative testing

% this would be done in a error basis level, how much error does the entire
% pose association have, in percentage, and in histogram bars plotting the
% data.

quantitativeTestData = quantitativeTest(gtDenseAssc, denseAssc);
% -- print the data
noAsscGt = length(find(quantitativeTestData == 2));
noAsscAUT = length(find(quantitativeTestData == 3));
totalGtAssc = nnz(denseAssc);
less1percent = length(find(quantitativeTestData < 0.01));
less5percent = length(find(quantitativeTestData < 0.05));
less20percent = length(find(quantitativeTestData < 0.2));

fprintf(1,'Poses that dont have association from GT: %d\n',noAsscGt);
fprintf(1,'Poses that dont have association from AUT:%d\n', noAsscAUT);
fprintf(1,'Associations that have less than 1 percent error: %d percent\n',less1percent);
fprintf(1,'Associations that have less than 5 percent error: %d percent\n', less5percent);
fprintf(1,'Associations that have less than 20 percent error: %d percent\n',less20percent);


%% visualize data
% g2oOpFile1 = 'input/db1-optimised.g2o';
% g2oOpFile2 = 'input/db2-optimised.g2o';
lA1 = landmarkAssc > 0;
lA2tmp = landmarkAssc(landmarkAssc>0);
lA2 = zeros(size(landmarks2,2));
lA2(lA2tmp) = 1;

% [vertices1,vCount1,landmarks1, lCount1, edges1, eCount1, ledges1, leCount1, dim1] =  readLandmarkG2oFile(g2oOpFile1, 'fc.txt', 1);
% [vertices2,vCount2,landmarks2, lCount2, edges2, eCount2, ledges2, leCount2, dim2] =  readLandmarkG2oFile(g2oOpFile2, 'fc.txt', 1);
depth_z = 50;
x_offset = 500;
y_offset = 0;
color1 = jet(size(vertices1,2));
plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, color1,'red');
color2 = getColor2(color1, size(vertices2,2), denseAssc);
% plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, 'blue','red');
plotGraph(vertices2, landmarks2, lA2, x_offset, y_offset, depth_z, color2,'cyan');
plotDenseMap(denseAssc, vertices1, vertices2, landmarks1, landmarks2, x_offset,y_offset, depth_z);

%% export both graphs to a single g2o
outFileName = 'outFile.g2o';

% -- sparsify ground truth
gtSparseAssc = gtDenseAssc;
totalPoses = length(gtSparseAssc);
i = 1;
infoCount = 0;
while(i<=totalPoses)
    if(gtSparseAssc(i) ~= 0)
        infoCount = infoCount + 1;
        if (mod(infoCount,5) ~= 0)
            gtSparseAssc(i) = 0;
        end
    end
    i = i + 1;
end

% exportTwinG2oFile(outFileName, vertices1, landmarks1, edges1, ledges1, vertices2, landmarks2, edges2, ledges2, gtDenseAssc);
% exportTwinG2oFile(outFileName, vertices1, landmarks1, edges1, ledges1, vertices2, landmarks2, edges2, ledges2, gtSparseAssc);
exportTwinG2oFile(outFileName, vertices1, landmarks1, edges1, ledges1, vertices2, landmarks2, edges2, ledges2, denseAssc);

%{
A = makeIncidenceMatrix(vertices,vCount,landmarks, lCount, edges,eCount, ledges, leCount);
distn1 = heatEmbeddingDistance(A, t_scale, vertices, landmarks, edges, ledges);
vertices1 = vertices;
landmarks1 = landmarks;
edges1 = edges;
ledges1 = ledges;

load(file2);
A = makeIncidenceMatrix(vertices,vCount,landmarks, lCount, edges,eCount, ledges, leCount);
distn2 = heatEmbeddingDistance(A, t_scale, vertices, landmarks, edges, ledges);
vertices2 = vertices;
landmarks2 = landmarks;
edges2 = edges;
ledges2 = ledges;


landmarkAssc = findLandmarkAssociations(file1, file2);
toc;
% start the seed growing algorithm. You have the both the descriptors, as
% nTscales x (nvertices x nvertices) so we consider the (i,j)th vertex, which
% gives us the heat at vertex i due to vertex j. and we consider the
% scales (here we have 5 of them).

% now that descriptors are defined, we can start off with the known
% associations. We need to keep a list of vertices, which we know are
% connected. We need to have a list of vertices, which we know are already
% associated with a certain vertex of the other graph.

% after we have the descriptor for that vertex
denseAssc = getDenseAssociations(distn1, distn2, landmarkAssc, vertices1, vertices2, landmarks1, landmarks2, edges1, edges2, ledges1, ledges2);

%}


