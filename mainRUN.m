clear;
clc;


%% Initialize loop

% This is to run the kitti dataset 100 times, in each we have to:
% 1. Load the initial file (optimised)
% 2. Generate a set of landmarks on it.
% 3. Copy the landmarks to unoptimised, save both.
% 4. call my algorithm
% 5. call just ICP


%% init

% -- total loops
totalLoopCount = 5;

% -- unoptimised file
fileNOp = 'datasets/vp/unoptimised/complete/victoriaPark-full.g2o';
fileNOpMat = 'datasets/vp/unoptimised/complete/victoriaPark-full.mat';

% -- optimised file
fileOpComplete = 'datasets/vp/optimised/complete/victoriaPark-fullOp.g2o';
fileOpCompleteMat = 'datasets/vp/optimised/complete/victoriaPark-fullOp.mat';

% -- ateFile AUT
fileAutResults = 'vpAUTResults.txt';
fileICPResults = 'vpICPResults.txt';

% -- path for results
fileNameOwlBaseWithPath = 'datasets/vp/optimised/split/victoriaPark-op-Dataset';
fileNameNOwlBaseWithPath = 'datasets/vp/unoptimised/split/victoriaPark-nop-Dataset';
joinedFilesDir = 'datasets/VPoutput/joined/';

% -- path for quantitative result
ateAUTFile = 'AUTateFile.txt';
ateICPFile = 'ICPateFile.txt';

%% loop through the files
% noiseLevels = zeros(1,totalLoopCount);

for i = 1:totalLoopCount
    noiseLevels(i) = i/10;
end
%  for i = 1:totalLoopCount
%      % -- fileName
%      joinedFileName = sprintf('victoriaPark-%d',i);
%      joinedFileNameWithPath = [joinedFilesDir,joinedFileName];
%      % -- generate noise levels
%      noiseLevels(i) = i/10;
%  end

% we just consider one noiseLevel.
% noiseLevels = 0.4;
joinedFileName = sprintf('victoriaPark-vNoise-');
joinedFileNameWithPath = [joinedFilesDir,joinedFileName];
AUTMultirunDriver(fileNameNOwlBaseWithPath, fileNameOwlBaseWithPath, joinedFileNameWithPath, ateAUTFile, ateICPFile, noiseLevels);
%     % -- print iteration
%     fprintf('Iteration: %d out of %d\n',i,totalLoopCount);
%     % -- run the code to run our algorithm on it
%     joinedFileName = sprintf('victoriaPark-%d',i);
%     joinedFileNameWithPath = [joinedFilesDir,joinedFileName];
   %AUTrun(fileNameNOwlBaseWithPath, fileNameOwlBaseWithPath, joinedFileNameWithPath, ateAUTFile, ateICPFile, i/10);
%     AUTMultirunDriver(fileNameNOwlBaseWithPath, fileNameOwlBaseWithPath, joinedFileNameWithPath, ateAUTFile, ateICPFile, i/10);
%  end
