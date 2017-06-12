function AUTMultirunDriver(unoptimisedFileBaseName, optimisedFileBaseName, outFileNameBase, atePrecisionFile, atePrecisionFileICP, vertexNoiseLevel)


file1 = [unoptimisedFileBaseName,'-1.mat'];
file2 = [unoptimisedFileBaseName,'-2.mat'];

landmarkAssc = findLandmarkAssociations(file1, file2);
[gtPoseAssc, gtDenseAssc] = findPoseAssociations(file1, file2);

% WE DO NOT ADD NOISE TO ASSOCIATIONS, WE JUST CHECK FOR POSE NOISE

% -- add noise to landmarkAssc
%  nlac = nnz(landmarkAssc);

% % -- add 6 landmark Noise
% % we should test for 4 noise levels
% % percent: 5,10,15,20, each time we run 2 tests with different landmark
% % pose noise:
% ln = randi([1,nlac],3,2);
%
%
% % -- add the noise
% for i = 1:3
%     tmp              = landmarkAssc(ln(i,1));
%     landmarkAssc(ln(i,1))  = landmarkAssc(ln(i,2));
%     landmarkAssc(ln(i,2)) = tmp;
% end

% -- run the AUTMultiRun for multiple calls
%  ofnb = outFileNameBase;
%  for i = 1:4
%      changeIndicesCount = ceil((nlac*landmarkAsscNoise)/2);
%      ln = randi([1,nlac],changeIndicesCount,2);
%      % -- add the noise
%      for j = 1:changeIndicesCount
%          tmp              = landmarkAssc(ln(j,1));
%          landmarkAssc(ln(j,1))  = landmarkAssc(ln(j,2));
%          landmarkAssc(ln(j,2)) = tmp;
%      end

%      suffix = sprintf('-l%d',i);
%      outFileNameBase = [ofnb,suffix];
    AUTMultirun(unoptimisedFileBaseName, optimisedFileBaseName, outFileNameBase, ...
        atePrecisionFile, atePrecisionFileICP, vertexNoiseLevel, landmarkAssc, gtPoseAssc, gtDenseAssc);
%  end

end
