function [R, T] = icpForLandmarkAssc(landmarkAssc, landmarks1, landmarks2)

%% Convert landmarks to lData
lCount1 = size(landmarks1, 2);
lCount2 = size(landmarks2, 2);
lData1 = zeros(3,lCount1);
lData2 = zeros(3,lCount2);
for i = 1:lCount1
    lData1(1,i) = landmarks1(i).id;
    lData1(2,i) = landmarks1(i).x;
    lData1(3,i) = landmarks1(i).y;
end
for i = 1:lCount2
    lData2(1,i) = landmarks2(i).id;
    lData2(2,i) = landmarks2(i).x;
    lData2(3,i) = landmarks2(i).y;
end

%% convert landmarkAssc to matchesMade

%matchesMade is a 2 column array which has the information of the assc
lAsscCount = length(landmarkAssc);
mmCount = nnz(landmarkAssc);
matchesMade = zeros(mmCount,2);
idx = 1;
for i = 1:lAsscCount
    if(landmarkAssc(i) ~= 0)
        matchesMade(idx,1) = i;
        matchesMade(idx,2) = landmarkAssc(i);
        idx = idx + 1;
    end
end

%% icp for matchesMade

[R, T] = icpForMatchesMade(matchesMade, lData1, lData2);

end


    