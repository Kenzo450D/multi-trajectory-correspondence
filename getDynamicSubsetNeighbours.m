function [indexEnds] = getDynamicSubsetNeighbours(vMap, vData, lData, lAssc, landmarkNLimit)
%GETDYNAMICSUBSETNEIGHBOURS Returns a dynamic subset of neighbours,
%depending on how spare the landmarks are spread across on the map.
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

%% start each point at the landmark points

indexEnds = zeros(size(lAssc,2),2);
lcData = lData(1,lAssc);
vIdx = vData(1,:);
vCount = size(vData,2);
% -- initialize
for i = 1:length(lcData)
    diffVIdx = vIdx - lData(1,lAssc(i));
    tmpPosIdx = find(diffVIdx > 0);
    tmpNegIdx = find(diffVIdx < 0);
    if (isempty(tmpPosIdx))
        endIdx = vIdx(end);
    else
        endIdx = tmpPosIdx(1);
    end
    if (isempty(tmpNegIdx))
        startIdx = vIdx(1);
    else
        startIdx = tmpNegIdx(end);
    end
    indexEnds(i,1) = startIdx;
    indexEnds(i,2) = endIdx;
end
initIndexEnds = indexEnds;
% fprintf(1,'Start: \n');
% disp(indexEnds);

% -- increase seed either side

% the two flags endflag and startflag are used to that if landmarks are
% really close to one another, the territory of a landmark would still not
% spill over to the other landmark. And provide bizaare results. This is
% the solution for a bug which made the second last territory to be more
% than 2000 units wide, even though it wasn't supposed to happen.

changeFlag = 1;
while(changeFlag == 1)
    changeFlag = 0;
    for i = 1:length(lcData)
        endFlag = 0;
        startFlag = 0;
        tmpEndIdx = indexEnds(i,2);
        findIdx = find(indexEnds(:,1) == tmpEndIdx);
        if (~isempty(findIdx) && i~=length(lcData))
            endFlag = 1;
        end
        tmpStartIdx = indexEnds(i,1);
        findIdx = find(indexEnds(:,2) == tmpStartIdx);
        if (~isempty(findIdx) && i~=1)
            startFlag = 1;
        end
        tmpEndIdx = indexEnds(i,2) + 1;
        findIdx = find(indexEnds(:,1) == tmpEndIdx);
        if (isempty(findIdx) && tmpEndIdx <= vCount && endFlag ~= 1)
            changeFlag = 1;
            indexEnds(i,2) = indexEnds(i,2) + 1;
        end
        tmpStartIdx = indexEnds(i,1) - 1;
        findIdx = find(indexEnds(:,2) == tmpStartIdx);
        if (isempty(findIdx) && tmpStartIdx > 0 && startFlag ~= 1)
            changeFlag = 1;
            indexEnds(i,1) = indexEnds(i,1) - 1;
        end 
    end
end

%% handle the first and last indexes

% There are lot of unneded associations which get made if not handled, and
% most of them are incorrect, it is wise to correct them. We choose a
% neighbourhood of 50 when such condition arises
maxNeighborhoodBorder = 50;

if (indexEnds(length(indexEnds),2) > initIndexEnds(length(indexEnds),2) + maxNeighborhoodBorder)
    indexEnds(length(indexEnds),2) = initIndexEnds(length(indexEnds),2) + maxNeighborhoodBorder;
end

if (indexEnds(1,1) < (initIndexEnds(1,1) - maxNeighborhoodBorder))
    indexEnds(1,1) = initIndexEnds(1,1) - maxNeighborhoodBorder;
end

% fprintf(1,'End: \n');
% disp(indexEnds);

%% handle small territories, by considering the information from immediate neighbourhood of the landmarks

% -- not required.

end