function makeKittiDataset(inputFileName, outputFileBase)

%% params
%fileName = 'input/outKitti.g2o';
%outputFileBase = 'kitti00Dataset';
%% victoria park params
%fileName = 'input/optimised/victoriaPark-fullOp.g2o';
%outputFileBase = 'victoriaParkOpDataset';

%% read file
[vertices,vCount,landmarks, lCount, edges,eCount, ledges, leCount, dim] =  readLandmarkG2oFile(inputFileName , 'fc.txt', 1);

%% call function
makeDataset( vertices, landmarks, edges, ledges, outputFileBase );