leData = zeros(7,sum(keepLandmark));
leCount = 0;
vData = vData2;
% for i = 1:eCount
%     eData(1,i) = vData(1,eData(1,i));
%     eData(2,i) = vData(1,eData(2,i));
% end
%% get the edges for all the vertices which are in contact with the landmark

% we should get the vertex id's for the landmark indexes. To do that, we find
% out, exactly where the landmark gets used first, and we declare the index
% right before it happens.

landmarkIdx = zeros(lCount, 1);
vData2 = vData;  % this is the vData with old vertex information
lCount = length(randomPoints);
lIdx = 0;

for i = 1:vCount
    for j = 1:lCount
        lpx = vData(2,i) - randomPoints(j,1);
        lpy = vData(3,i) - randomPoints(j,2);
        lpDist = sqrt(lpx*lpx + lpy*lpy);
        if (lpDist < landmarkLookupRange * avgPseStep)
            % -- fill up the edge dx,dy
            leCount = leCount + 1;
            leData(1,leCount) = i;
            leData(2,leCount) = j;
            leData(3,leCount) = lpx;
            leData(4,leCount) = lpy;
            % -- fill up the information matrix
            leData(5,leCount) = 1.58114;
            leData(7,leCount) = 1.58114;
            % -- fill up the landmarkIdx
            if (landmarkIdx(j) == 0)
                tmpVarIDx = vData(1,i) + 1;
                % -- find the index in the landmarks
                foundChk = find(landmarkIdx == tmpVarIDx);
                while(~isempty(foundChk))
                    tmpVarIDx = tmpVarIDx + 1;
                    foundChk  = find(landmarkIdx == tmpVarIDx);
                end
                landmarkIdx(j) = tmpVarIDx;
                % -- update vData indexes
                for k = (i+1):vCount
                    vData(1,k) = vData(1,k)+1;
                end
                % -- add information to lData
                lIdx = lIdx + 1;
                lData(1,lIdx) = landmarkIdx(j);
                lData(2,lIdx) = randomPoints(j,1);
                lData(3,lIdx) = randomPoints(j,2);
            end
        end
    end
end