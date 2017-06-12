function [indexEnds] = getSubsetNeighbours(vMap, vData, lData, lAssc, landmarkNLimit)
%GETSUBSETNEIGHBOURS Returns a subset of neighbours of landmarks, for each landmark.
% Input:
%  vMap: inverse vertex map to get original vertex indices
%  vData: matrix form of vertex information
%  lData: matrix form of landmark information
%  lAssc: landmark Association indices
%  landmarkNLimit: number of landmarks to be chosen
% Output:
%  indexEnds: index based on vData, of limits on subsets, just has the start
%             and end index

% Note: We would only use the indices in lAssc to form the indexEnds

%% loop over lAssc, find out each limit
limitEitherSide = landmarkNLimit/2;
vIdx = vData(1,:);
indexEnds = zeros(size(lAssc,2),2);
for i = 1:length(lAssc)
    diffVIdx = vIdx - lData(1,lAssc(i));
    tmpPosIdx = find(diffVIdx > 0);
    tmpNegIdx = find(diffVIdx < 0);
    if length(tmpPosIdx >= limitEitherSide)
        endIdx = tmpPosIdx(limitEitherSide);
    else
        endIdx = tmpPosIdx(end);
    end
    if length(tmpNegIdx >= limitEitherSide)
        startIdx = tmpNegIdx(length(tmpNegIdx) - limitEitherSide);
    else
        startIdx = tmpNegIdx(1);
    end
    indexEnds(i,1) = startIdx;
    indexEnds(i,2) = endIdx;
end
end