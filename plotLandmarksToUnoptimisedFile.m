function plotLandmarksToUnoptimisedFile(fileNameOptimised, fileNameUnoptimised, fileNameLandmark, outFileName)
%% initialization
%fileNameUnoptimised = 'input/unoptimised/kitti00.g2o';
%fileNameOptimised = 'input/optimised/kitti00op.g2o';
%fileNameLandmark = 'input/optimisedLandmark/kitti00opWl.g2o';
filesCreatedFile = 'fc.txt';
%outFileName = 'outKitti.g2o';

%% read the files
% ---- read the unoptimised file (the file to add landmarks to)
[vertices1, vCount1, landmarks1, lCount1, edges1,eCount1, ledges1, leCount1, ~] =  readLandmarkG2oFile(fileNameUnoptimised, filesCreatedFile, 1);
edges1 = removeLoopClosures(edges1);
eCount1 = size(edges1,2);
% ---- load the optimised file
[verticesOptimised, vCountOptimised, landmarksOptimised, lCountOptimised, edgesOptimised, eCountOptimised, ledgesOptimised, leCountOptimised, ~] =  readLandmarkG2oFile(fileNameOptimised, filesCreatedFile, 1);
edgesOptimised = removeLoopClosures(edgesOptimised);
eCountOptimised = size(edgesOptimised,2);
% ---- load the optimised file with the landmarks
[verticesOpwl, vCountOpwl, landmarksOpwl, lCountOpwl, edgesOpwl, eCountOpwl, ledgesOpwl, leCountOpwl, dim3] =  readLandmarkG2oFile(fileNameLandmark, filesCreatedFile, 1);

% -- convert to matrices
[vData1, lData1, eData1, leData1] = getMatrixForm(vertices1, landmarks1, edges1, ledges1);
[vDataOp, lDataOp, eDataOp, leDataOp] = getMatrixForm(verticesOptimised, landmarksOptimised, edgesOptimised, ledgesOptimised);
[vDataOpWl, lDataOpWl, eDataOpWl, leDataOpWl] = getMatrixForm(verticesOpwl, landmarksOpwl, edgesOpwl, ledgesOpwl);


%% create a map of the pose indexes
% ---- find total number of pose indexes along with landmarks
totalPoseCountOpwl = vCountOpwl + lCountOpwl;
poseMap  = zeros(totalPoseCountOpwl,1);
vertexMap = zeros(vCountOpwl);

% ---- map the pose index of "optimised" file and "optimised with landmark" file
for i = 1:vCountOpwl
    vertexMap(i) = vDataOpWl(1,i);
end

%% change the pose indexes of the unoptimised files
% -- take a backup
vDataInit = vData1;
lDataInit = lData1;
eDataInit = eData1;
leDataInit = leData1;

% -- convert the vData indexes
for i = 1:vCount1
    vData1(1,i) = vertexMap(vData1(1,i));
end

% -- convert the edge indexes
for i = 1:eCount1
    eData1(1,i) = vertexMap(eData1(1,i));
    eData1(2,i) = vertexMap(eData1(2,i));
end

%% add the landmarks and export the file
leData1 = leDataOpWl;
lData1 = lDataOpWl;


exportLandmarkG2oFile(outFileName, vData1, lData1, eData1, leData1);

end

