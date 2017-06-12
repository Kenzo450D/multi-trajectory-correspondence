function [denseAssc] = getDenseAssociations(distn1, distn2, landmarkAssc, vertices1, vertices2, landmarks1, landmarks2, edges1, edges2, ledges1, ledges2)
%GETDENSEASSOCIATIONS Get a list of dense associations given a list of
%landmark Associations, and the vertex and edge informations.
% Input:
%   distn1 : heat distance descriptors for graph1
%   distn2 : heat distance descriptros for graph2
%   landmarkAssc: landmark associations, a binary array of size of lCount1
%   vertices1: vertices set of graph1
%   vertices2: vertices set of graph2
%   landmark1: landmarks set of graph1
%   landmark2: landmarks set of graph2
%   edges1: edges set of graph1
%   edges2: edges set of graph2
%   ledges1: landmark edges set of graph1
%   ledges2: landmark edges set of graph2

%% initializations
vCount1 = size(vertices1, 2);
vCount2 = size(vertices2, 2);
lCount1 = size(landmarks1, 2);
lCount2 = size(landmarks2, 2);
totalVertexCount1 = vCount1 + lCount1;
totalVertexCount2 = vCount2 + lCount2;

% - debug
fprintf(1,'In dense Association function: ');
fprintf(1,'Graph1: ');
fprintf(1,'    vCount: %d',vCount1);
fprintf(1,'    lCount: %d', lCount1);
fprintf(1,'Graph2: ');
fprintf(1,'    vCount: %d',vCount2);
fprintf(1,'    lCount: %d', lCount2);

%% convert vertex and edge data to matrix
[vData1, lData1, eData1, leData1] = getMatrixForm(vertices1, landmarks1, edges1, ledges1);
[vData2, lData2, eData2, leData2] = getMatrixForm(vertices2, landmarks2, edges2, ledges2);

%% initialize the dense association matrix
denseAssc = zeros(totalVertexCount1,1);
alreadyAssc1 = zeros(totalVertexCount1,1);
alreadyAssc2 = zeros(totalVertexCount2,1);

% -- start off with the first association

landmarkIdx = 1;
index2 = 1;

while(true)
    % -- take first landmark correspondence
    if (landmarkIdx > length(landmarkAssc))
        break;
    end
    if (landmarkAssc(landmarkIdx) == 0)
        landmarkIdx = landmarkIdx + 1;
        continue;
    end
    % -- get it's neighbours
    [vertexNeighbours1] = getNeighbours(leData1, lData1,landmarkIdx);
    [vertexNeighbours2] = getNeighbours(leData2, lData2,landmarkAssc(landmarkIdx));

    % -- check for empty neighbour list
    if (isempty(vertexNeighbours1) || isempty(vertexNeighbours2))
        fprintf(1,'Empty neighbour list either for: ');
        fprintf(1,'Graph1: %d\tGraph2: %d\n',lData2(1,landmarkIdx), lData2(1,landmarkAssc(landmarkIdx)));
        landmarkIdx = landmarkIdx + 1;
        continue;
    end
    
    % -- check if the neighbours haven't already been associated
    % ** -- TODO here -- **
    alreadyAsscIdx = denseAssc(vertexNeighbours1);
    delFlag1 = [];
    delFlag2 = [];
    for i = 1:length(alreadyAsscIdx)
        if (alreadyAsscIdx(i) ~= 0)
            delFlag1 = [delFlag1,i];
            delFlag2 = [delFlag2,alreadyAsscIdx(i)];
        end
    end
    if (~isempty(delFlag1))
        fprintf(1,'landmarkIdx: %d\n',landmarkIdx);
        vertexNeighbours1(delFlag1) = [];
        vertexNeighbours2(ismember(vertexNeighbours2, delFlag2)) = [];
    end
    
    % -- take a backup of unasscoaited neighbours (if required later)
    initVertexNeighbours1 = vertexNeighbours1;
    initVertexNeighbours2 = vertexNeighbours2;
    vertexIdx1_length = length(vertexNeighbours1);
    vertexIdx2_length = length(vertexNeighbours2);
    
    % -- print the details for debug
    fprintf(1,'Graph1 Landmark: %d\n',lData1(landmarkIdx));
    fprintf(1,'\t\tNeighbours: ');
    for i = 1:length(vertexNeighbours1)
        fprintf(1,'%d ',vertexNeighbours1(i));
    end
    fprintf(1,'\n');
    fprintf(1,'Graph2 Landmark: %d\n',lData2(landmarkAssc(landmarkIdx)));
    fprintf(1,'\t\tNeighbours: ');
    for i = 1:length(vertexNeighbours2)
        fprintf(1,'%d ',vertexNeighbours2(i));
    end
    fprintf(1,'\n');
    % -- end debug
    
    
    % -- loop through all neighbour combinations to find out the best
    % matches & save matches in the dense association matrix
    for i = 1:vertexIdx1_length
        desc1 = squeeze(distn1(:,vertexNeighbours1(i),(vCount1+landmarkIdx)));
        if (~isempty(vertexNeighbours2))
            desc2 = squeeze(distn2(:,vertexNeighbours2(1),(vCount2+landmarkAssc(landmarkIdx))));
            minVal = desc1 - desc2;
            minValNorm = norm(minVal);
            vertexIdx2_length = length(vertexNeighbours2);
            bestMatch = vertexNeighbours2(1);
            bestMatchIdx = 1;
            for j = 2:vertexIdx2_length;
                desc2 = squeeze(distn2(:,vertexNeighbours2(j),(vCount2+landmarkAssc(landmarkIdx))));
                tmpMinVal = desc1 - desc2;
                tmpMinValNorm = norm(tmpMinVal);
                if (tmpMinValNorm < minValNorm)
                    minValNorm = tmpMinValNorm;
                    bestMatch = vertexNeighbours2(j);
                    bestMatchIdx = j;
                end
            end
            % -- save the best match
            denseAssc(vertexNeighbours1(i)) = bestMatch;
            fprintf('Matched: Graph1: %d\twith\tGraph2: %d\n',vertexNeighbours1(i), bestMatch);
            % -- delete the element from the second set of neighbours
            vertexNeighbours2(bestMatchIdx) = [];
        end
    end
    % -- go to next landmark
    landmarkIdx = landmarkIdx + 1;
end

end