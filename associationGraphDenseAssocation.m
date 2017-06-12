function [matchesMade] = associationGraphDenseAssocation(agVertexList, eucDiff1, eucDiff2, totalVertexCount1, vertices1, vertices2, vd1, vd2, landmarks1, landmarks2, landmarkAssc)

%% init the threshold
hkdThresh = 2;

%% compute the edge constraints of the association graph

agVertexCount = size(agVertexList,1);

tic;
agEdgeCount = 0;
agList = zeros(agVertexCount,agVertexCount); %adjacency List
agAM = zeros(agVertexCount,agVertexCount); %adjacency matrix
%agLogical = logical(agVertexCount,
agListEndIndexes = ones(agVertexCount,1);
for i = 1:agVertexCount
    v1_1Idx = agVertexList(i,1); % vertex of graph1, in vertex1
    v1_2Idx = agVertexList(i,2); % vertex of graph2, in vertex1
    for j = 1:agVertexCount
        if (i ~= j)
            v2_1Idx = agVertexList(j,1); % vertex of graph1, in vertex2
            v2_2Idx = agVertexList(j,2); % vertex of graph2, in vertex2
            distDiff1 = abs(eucDiff1(v1_1Idx, v2_1Idx) - eucDiff2(v1_2Idx,v2_2Idx));
            %distDiff2 = getDistanceDiffVertex(vData1, vData2, v1_1Idx, v1_2Idx,v2_1Idx, v2_2Idx);
            agAM(i,j) = distDiff1;
            if (distDiff1 < 2)
                agAM(i,j) = totalVertexCount1;
                agList(i,agListEndIndexes(i)) = j;
                agListEndIndexes(i) = agListEndIndexes(i) + 1;
                agEdgeCount = agEdgeCount + 1;
            end
        end
    end
end
fprintf(1,'Time taken to create edges: %d\n',toc);
fprintf(1,'Number of vertices: %d\n',agVertexCount);
fprintf(1,'Total Edges made: %d\n', agEdgeCount);

%% find the largest clique

meanVal = mean(mean(agAM));
%weightsMap = exp(-(1.5*agAM)/meanVal);
weightsMap = exp(-(1.3*agAM)/meanVal);
weightsMap(1:agVertexCount+1:agVertexCount*agVertexCount) = 0;

% -- SVD of the matrix
[v,d] = eig(weightsMap);

%{
% Notes:
% For each candidate assignment a = (i,i') there is an associated score or
% affinity that measures how well feature i belongs to P, i' belongs to Q.
% We have a constraint to allow one feature from P to match at most one
% feature from Q. Given a list L of n candidate assignments we store the
% affinities on every assignment a belongs to L, and every pair of
% assignments a,b belongs to L in the n x n matrix M (here weightsMap)

% L here is agVertexList

% 2. Let x* be the principal eigenvector of M. Initialize the solution
% vector x with the nx1 zero vector. Initialize L with the set of all
% candidiate assignments.

% 3. Find a* = argmax a∈L (x*(a)). If x*(a*) = 0,  stop and return the
% solution.
%}

%% consider the first eigenvector of weightsMap

% we consider the first column of V, which would correspond to the first
% eigen vector.

%% find eigenvector corresponding to the largest eigenvalue
[~,largestEigValIdx] = max(diag(d));
eigVector = v(:,largestEigValIdx);

% consider the absolute values of the eigenvector
eigVector = abs(eigVector);

% find the end-counter for the values (See if Required)
medVal = median(eigVector);
tmp = abs(eigVector-medVal);
[~,medIdx] = min(tmp);

% sort the eigenVectors and consider only the maximum ones
maxMatchCounter = agVertexCount/2;
currentMatchCounter = 0;
matchesMade = [];
i = 1;
% -- Get threshold
[minVal,minIdx] = min(eigVector);
%eigVector(minIdx) = medVal;
%[minVal2] = min(eigVector);
thresholdStopVal = (medVal + minVal)/2;
thresholdStopVal = medVal;
%eigVector(minIdx) = minVal;



% -- loop starts here
while(currentMatchCounter < maxMatchCounter || i <= agVertexCount)
    [maxVal,maxIdx] = max(eigVector);
    % ---- check if value is less than median of the eigenvector
    if (maxVal <= thresholdStopVal)
        break;
    end
    % ---- get selected Pair
    selPair_idx1 = agVertexList(maxIdx,1);
    selPair_idx2 = agVertexList(maxIdx,2);
    % check for heat diffusion comparison less than threshold
    cmpHd = vd1(selPair_idx1) - vd2(selPair_idx2);
    if (norm(cmpHd) > hkdThresh)
        continue;
    end
    % check for collisions for selected value
    if (isempty(matchesMade))
        tmp1 = [];
        tmp2 = [];
    else
        tmp1  = find(matchesMade(:,1) == selPair_idx1);
        tmp2 = find(matchesMade(:,2) == selPair_idx2);
    end
    if (isempty(tmp2) && isempty(tmp1))
        % ---- pair doesn't collide, add the pair
        tmpMatch = [selPair_idx1, selPair_idx2];
        matchesMade = [matchesMade;tmpMatch];
        currentMatchCounter = currentMatchCounter + 1;
    end
    % ---- make current max value as min value
    eigVector(maxIdx) = minVal;
    % ---- display the results intermittently
    % display the results every few iterations
    if ( mod(i,1000) == 0)
        %denseAssc = zeros(totalVertexCount1,1);
        %if (~isempty(matchesMade))
            %denseAssc(matchesMade(:,1)) = matchesMade(:,2);
            %displayResults(denseAssc, landmarkAssc, landmarks1, landmarks2, vertices1, vertices2);
        %end
        %hold off; figure();
        fprintf(1,'i = %d\n',i);
    end
    % ---- next loop iterator
    i = i+1;
end


% denseAssc = zeros(totalVertexCount1,1);
% if (~isempty(matchesMade))
%     denseAssc(matchesMade(:,1)) = matchesMade(:,2);
% end



end