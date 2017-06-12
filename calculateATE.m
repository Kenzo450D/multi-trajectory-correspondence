function [ate] = calculateATE(vData1, vData2, matchesMade)
%CALCULATEATE Calculates ATE in the matches made
% -- 


% To calculate the ATE, we compare the values of the gtMatch indexes to
% that of the aut Match indexes. 

% -- Convert gtPoseAssc to matchesMade

sumDiffn = 0;
for i = 1:length(matchesMade)
    v1x = vData1(2,matchesMade(i,1));
    v1y = vData1(3,matchesMade(i,1));
    v2x = vData2(2,matchesMade(i,2));
    v2y = vData2(3,matchesMade(i,2));
    % -- euclidean difference
    diffn = sqrt((v1x - v2x)^2 + (v1y - v2y)^2);
    % -- add the diff
    sumDiffn = sumDiffn + diffn;
end

ate = sumDiffn / length(matchesMade);
end