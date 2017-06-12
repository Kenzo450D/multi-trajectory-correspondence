function AUTMultirun(unoptimisedFileBaseName, optimisedFileBaseName, outFileNameBase, atePrecisionFile, atePrecisionFileICP, vertexNoiseLevel, landmarkAssc, gtPoseAssc, gtDenseAssc)
%AUTMULTIRUN Runs the same set of input multiple times with different
%initializations of noise in pose and landmark noise. 2 times each, thereby
%resulting in a 4 times iteration.
% landmarkAssc: provide noisy landmarkAssc
% gtPoseAssc:  gtDenseAssc, vertices1, vertices2)
%% initialize
file1 = [unoptimisedFileBaseName,'-1.mat'];
file2 = [unoptimisedFileBaseName,'-2.mat'];
optimisedFile1 = [optimisedFileBaseName,'-1.mat'];
optimisedFile2 = [optimisedFileBaseName,'-2.mat'];

%% get sparse landmark association

% We consider this as input

%{
% -- Victoria park
% file1 = 'input/databaseFile-1.mat';
% file2 = 'input/databaseFile-2.mat';
% optimisedFile1 = 'input/optimised/victoriaParkOpDataset-1.mat';
% optimisedFile2 = 'input/optimised/victoriaParkOpDataset-2.mat';

% -- kitti dataset
%file1 = 'input/unoptimised/kitti00Dataset-1.mat';
%file2 = 'input/unoptimised/kitti00Dataset-2.mat';


% landmarkAssc = findLandmarkAssociations(file1, file2);
% [gtPoseAssc, gtDenseAssc] = findPoseAssociations(file1, file2);
%}

%% t_scale parameters

%t_scale = [0.001,0.01,1,10,50];
%t_scale_indicator = [0.99,0.95,0.92,0.9,0.85, 0.8, 0.7, 0.6];
%t_scale_indicator = [0.99,0.95,0.92,0.9,0.85, 0.8];
%t_scale_indicator = [0.99,0.95,0.9,0.85, 0.8, 0.75, 0.7, 0.65,0.6];
t_scale_indicator = [0.99,0.95,0.9,0.85, 0.8, 0.75];%
n_perc_eigval_indx = 5/100;%0.5/100;
%t_scale = [0.001,0.01,1,10,100];


tic;
%% load each file to find
% =========================================================================
% FILE 1
% =========================================================================
% -- load the file
load(file1);
% -- save the non-noisy vertices
initV1 = vertices;

% -- add noise to each of the files
%[vertices] = addNoiseToTrajectoryPoses(vertices, vertexNoiseLevel);
%[vertices, avgOdomLength] = addSuccessiveNoise(vertices,vertexNoiseLevel);
%[landmarks] = addNoiseToLandmarkPoses(landmarks, 30, 2);

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
[distn1, vd1, tScale1] = heatEmbeddingDescriptor(L, t_scale_indicator, [], n_perc_eigval_indx, vertices, landmarks, edges, ledges, l1Assc);
vertices1 = vertices;
landmarks1 = landmarks;
edges1 = edges;
ledges1 = ledges;

% =========================================================================
% FILE 2
% =========================================================================
load(file2);
% -- save the non-noisy vertices
initV2 = vertices;

% -- add noise to each of the files
%[vertices] = addNoiseToTrajectoryPoses(vertices, vertexNoiseLevel);
%[vertices, avgOdomLength] = addSuccessiveNoise(vertices,vertexNoiseLevel);
%[landmarks] = addNoiseToLandmarkPoses(landmarks, 10, 2);

% ---- make the adjacency matrix
A = makeAdjacencyMatrix(vertices, landmarks, edges, ledges);
% ---- get the laplacian matrix
L = getlaplacianMatrix(A);

% -- prepare associations
% find the poses in map2 which have an association in map1
l2Assc = landmarkAssc>0;
l2Assc = landmarkAssc(l2Assc); % TODO: Check this

[distn2,vd2, tScale2] = heatEmbeddingDescriptor(L, t_scale_indicator, [], n_perc_eigval_indx, vertices, landmarks, edges, ledges, l2Assc);
vertices2 = vertices;
landmarks2 = landmarks;
edges2 = edges;
ledges2 = ledges;
toc;

[matchesMade, denseAssc] = getDenseAssociationsLocal(distn1, distn2, vd1, vd2, ...
    landmarkAssc, vertices1, vertices2,...
    landmarks1, landmarks2, ...
    edges1, edges2, ...
    ledges1, ledges2, gtPoseAssc, ...
    tScale1, tScale2);

landmarkNoiseLevelArray = [2,5,7,10,12,15,17,20];
repeatForTimes = 6;

load(optimisedFile1);
opVertices1 = vertices;
load(optimisedFile2);
opVertices2 = vertices;
initVData1 = getVertexMatrixForm(opVertices1);
initVData2 = getVertexMatrixForm(opVertices2);

% we require the work done to have the same vertex noise Level, while we
% change through landmark noise. Hence we initialize the vertex error here
% outside the loop and we do not change this while we are inside it.

for noiseIdx = 1:length(vertexNoiseLevel)
    load(file1);
    vertices1 = vertices;
    [vertices1] = addSuccessiveNoise(vertices1,vertexNoiseLevel(noiseIdx));
    load(file2);
    vertices2 = vertices;
    [vertices2] = addSuccessiveNoise(vertices2,vertexNoiseLevel(noiseIdx));
    
    % now we have different landmark noise levels, we would run each of them 5
    % times, get a an even distribution.
    for loopIdx = 1:5
        % -- we increase landmark Noise every iteration
        landmarkNoiseLevel = landmarkNoiseLevelArray(loopIdx);
        for repeatLoop = 1:repeatForTimes
            %% set up outFileName
            outFileSuffix = sprintf('%d-output',noiseIdx);
            outFileNameNewBase = [outFileNameBase,outFileSuffix];
            
            %% add noise to the trajectory
            load(file1);
            landmarks1 = landmarks;
            [landmarks1] = addNoiseToLandmarkPoses(landmarks1, landmarkNoiseLevel, 1);
            load(file2);
            landmarks2 = landmarks;
            [landmarks2] = addNoiseToLandmarkPoses(landmarks2, landmarkNoiseLevel, 1);
            
           %% run Icp to the dataset
            [Ricp1, Ticp1, ER1, t1, newMatchesMade, denseAssc]  = runICP(matchesMade, vertices1, vertices2, landmarks1, landmarks2, vd1, vd2);
            
            %% plot
            %     cropPointT1v1 = min(newMatchesMade(:,1));
            %     cropPointT1v2 = max(newMatchesMade(:,1));
            %     cropPointT2v1 = min(newMatchesMade(:,2));
            %     cropPointT2v2 = max(newMatchesMade(:,2));
            %     opFileName1 = optimisedFile1;
            %     opFileName2 = optimisedFile2;
            %     nameFig = [outFileNameNewBase,'AUT'];
            %     plotTheCroppedGraphs(opFileName1, opFileName2, newMatchesMade, cropPointT1v1, cropPointT1v2, cropPointT2v1, cropPointT2v2, nameFig);
            
           %% export both graphs to a single g2o
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
            
            outFileName = [outFileNameNewBase,'.g2o'];
            exportTwinG2oFile(outFileName, vertices1, landmarks1, edges1, ledges1, vertices2, landmarks2, edges2, ledges2, denseAssc);
            
           %% Run code to get the ATE and Precision results on it.
            
            % -- convert gtPoseAssc to gtMatchesMade
            poseIdx = gtPoseAssc>0;
            posesToCompare = gtPoseAssc(poseIdx);
            gtMatchesMade = zeros(length(posesToCompare),2);
            for i = 1:length(posesToCompare)
                if (poseIdx(i) ~= 0)
                    gtMatchesMade(i,1) = i;
                    gtMatchesMade(i,2) = gtPoseAssc(i);
                end
            end
            
            vData1 = getVertexMatrixForm(vertices1);
            vData2 = getVertexMatrixForm(vertices2);
            
            [ate] = calculateATE(initVData1, initVData2, newMatchesMade);
            [p, relaxPrecision] = calculatePrecision(gtMatchesMade, newMatchesMade, (size(vData1,2)*0.001));
            
           %% save the data so that we can reload to plot again
            OutMatFilename = [outFileNameNewBase,'ForPlot.mat'];
            save(OutMatFilename,'denseAssc', 'vertices1', 'vertices2', 'edges1', 'edges2', 'landmarks1', 'landmarks2', 'ledges1', 'ledges2', 'landmarkAssc');
            
            %outFile = 'atePrecision.txt';
            fid = fopen(atePrecisionFile,'a')';
            
            fprintf(1,'AUT: %f %f %f\n',ate,p,relaxPrecision);
            fprintf(fid,'%f %f %f\n',ate,p,relaxPrecision);
            fclose(fid);
            
           %% ICP ONLY on the noisy information, to check the performance of ICP
            
            % we save the result on icp ate precision file
            outMatICPFileName = [outFileNameNewBase, '-ICPMatchesMade.mat'];
            %icpMainDataPassed(vertices1, vertices2, gtPoseAssc, atePrecisionFileICP);
            cropPointT1v1 = min(newMatchesMade(:,1));
            cropPointT1v2 = max(newMatchesMade(:,1));
            cropPointT2v1 = min(newMatchesMade(:,2));
            cropPointT2v2 = max(newMatchesMade(:,2));
            [icpMatchesMade] = icpSupervisedMainDataPassedWithCrop(vertices1, vertices2, ...
                landmarks1, landmarks2, cropPointT1v1, cropPointT1v2, ...
                cropPointT2v1, cropPointT2v2,   landmarkAssc, atePrecisionFileICP, ...
                outMatICPFileName);
            % [icpMatchesMade] = icpSupervisedMainDataPassed(vertices1, vertices2, landmarks1, landmarks2, gtPoseAssc, landmarkAssc, atePrecisionFileICP, outMatICPFileName);
            
            % -- print the results here
            %     cropPointT1v1 = min(icpMatchesMade(:,1));
            %     cropPointT1v2 = max(icpMatchesMade(:,1));
            %     cropPointT2v1 = min(icpMatchesMade(:,2));
            %     cropPointT2v2 = max(icpMatchesMade(:,2));
            %     nameFig = [outFileNameNewBase,'ICP'];
            %     plotTheCroppedGraphs(opFileName1, opFileName2, icpMatchesMade, cropPointT1v1, cropPointT1v2, cropPointT2v1, cropPointT2v2, nameFig);
            
            % -- take only those from icpMatchesMade which are present in newMatchesMade
            icpNMM = zeros(size(newMatchesMade));
            for i = 1:length(icpNMM)
                idx = newMatchesMade(i,2);
                idx2 = icpMatchesMade(:,2) == idx;
                icpNMM(i,:) = icpMatchesMade(idx2,:);
            end
            [ate] = calculateATE(initVData1, initVData2, icpNMM);
            
            % -- calculate precision
            [p, relaxPrecision] = calculatePrecision(gtMatchesMade, icpMatchesMade, (size(vData1,2)*0.001));
            
            outFile = atePrecisionFileICP;
            fid = fopen(outFile,'a')';
            
            fprintf(1,'ICP: %f %f %f\n',ate,p,relaxPrecision);
            fprintf(fid,'%f %f %f\n',ate,p,relaxPrecision);
            fclose(fid);
            
           %% save the matches made in the file
            save(outMatICPFileName,'icpMatchesMade');
            close all;
        end
    end
end
end