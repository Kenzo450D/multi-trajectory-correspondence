function [denseAssc] = getDenseAssociationsSubset(vd1, vd2, v1start, v1end, v2start, v2end, eucDiff1, eucDiff2, neighbourLimit, totalVertexCount1, tScale1, tScale2, poseNeighbourhoodMap1, poseNeighbourhoodMap2)
% -- Input
%   vd1 : vertex descriptor set 1
%   vd2 : vertex descriptor set 2
%   v1start : start boundary of vertex set 1
%   v1end : end boundary of vertex set 1
%   v2start : start boundary of vertex set 2
%   v2end : end boundary of vertex set 2
%   eucDiff1 : pose steps for each pose in graph 1
%   eucDiff2 : pose steps for each pose in graph 2
%   neighbourLimit: number of neighbours to be considered in Association Graph
%   totalVertexCount1: total vertex count for the first graph, needed to arrange data
%   tScale1 : tScale parameters for graph 1
%   tScale2 : tScale parameters for graph 2
%   *POSENEIGHBOURHOODMAP* : It is to denote how many neighbours does that
%   particular landmark territory have, it is so that a higher scale of
%   diffusion can be given higher importance in spare areas, and a lower
%   scale of diffusion can be given importance in dense areas.
%   poseNeighbourhoodMap1 : pose Neighbourhood information for graph 1
%   poseNeighbourhoodMap2 : pose Neighbourhood information for graph 2)
%   poseNeighbourhoodMap1

% -- count the number of poses

% Count the number of poses, if the pose count is small, we should give
% more importance to the low t_scale information of the heat kernel
% descriptors for the vertex comparisons. If the pose count is higher, then
% more information for the higher comparisons.

poseDiff = v1end - v1start;

% distnType contains information about the type of distance. 
% Type 0 to 5 would be based on the tScale steps.
    

% ---- make sparse matrix to make their adjacency graph

agVertexList = [];
%cmpHds = ones(usvCount1, usvCount2)*100;
% -- form the vertices
tic;
descr_threshold = 0.05;
for i = v1start:v1end
    index1 = i;
    hdv1   = vd1(:,:,index1);
    agVertexTmpList = [];
    agIndexTmpList = [];
    cmpHds_arr = [];
    %fprintf(1,'i = %d out of %d\n',i,usvCount1);
    for j = v2start:v2end
        index2 = j;
        hdv2   = vd2(:,:,index2);
        if (norm(hdv1) < 1e-6 || norm(hdv2) < 1e-6)
            continue;
        end
        % -- compare heat descriptors based on the difference of the
        % sections
        cmpHd = getComparedHeatDescriptor(hdv1, hdv2, tScale1, tScale2, poseNeighbourhoodMap1(i), poseNeighbourhoodMap2(j));
        %cmpHd = norm(hdv1 - hdv2); % compare between two heat descriptor for vertices
        if (cmpHd < descr_threshold)
            cmpHds_arr = [cmpHds_arr; cmpHd];
            tmpVal = [index1, index2];
            agVertexTmpList = [agVertexTmpList;tmpVal];
        end
    end
    [tval1 tindx1]=sort(cmpHds_arr);
     if length(cmpHds_arr)>neighbourLimit
        for j = 1:neighbourLimit
            agVertexList = [agVertexList; agVertexTmpList(tindx1(j),:)];
        end
     elseif length(cmpHds_arr)>0
            agVertexList = [agVertexList; agVertexTmpList];
     %else
     %       disp(sprintf('No Match Found for index: %d in left trajectory',index1));
     end
    %fprintf(1,'length of agVertexList: %d\n',length(agVertexList));
end

agVertexCount = length(agVertexList);
fprintf(1,'Time taken to create vertices: %d\n',toc);

if (agVertexCount == 0)
    denseAssc = zeros(totalVertexCount1,1);
    return;
end


%% only select the top 5 nearest neighbours (in terms of cmpHds) for each vertex

tic;
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
            if (distDiff1 < 0.001)
                agList(i,agListEndIndexes(i)) = j;
                agListEndIndexes(i) = agListEndIndexes(i) + 1;
            end
        end
    end
end
fprintf(1,'Time taken to create edges: %d\n',toc);
fprintf(1,'Number of vertices: %d\n',agVertexCount);
%% plot all the correspondences

%{
%  lA1 = landmarkAssc > 0;
%  lA2tmp = landmarkAssc(landmarkAssc>0);
%  lA2 = zeros(size(landmarks2,2));
%  lA2(lA2tmp) = 1;
%  depth_z = 50;
%  x_offset = 500;
%  y_offset = 0;
%  plotGraph(vertices1, landmarks1, lA1, edges1, ledges1, 0, 0, depth_z, 'blue','red');
%  plotGraph(vertices2, landmarks2, lA2, edges2, ledges2, x_offset, y_offset, depth_z, 'green','cyan');
%  plotDenseMap_twinIndex(agVertexList, vertices1, vertices2, landmarks1, landmarks2, x_offset,y_offset, depth_z);
%}

%% make a check of the largest connected component

% RESULT : there is just one connected component
%{
% Mark all the vertices as not visited
visitedV = zeros(agVertexCount,1);

% Create a queue for BFS
import java.util.LinkedList
qu = LinkedList();
%q.add('item1');
%q.add(2);
%q.add([3 3 3]);
%item = q.remove();
%q.add('item4');

% Mark the current node as visited and enqueue it



for i = 1:agVertexCount
    s = i;
    if (visitedV(s) == 1)
        continue;
    end
    visitedV(s) = 1;
    qu.add(s);
    while(qu.size() ~= 0)
        % Dequeue a vertex from queue and print it
        s = qu.remove();

        % Get all adjacent vertices of the dequeued vertex s
        % If a adjacent has not been visited, then mark it visited
        % and enqueue it

        for j = 1:(agListEndIndexes(s)-1)
            if(visitedV(agList(s,j)) == 0)
                visitedV(agList(s,j)) = 1;
                qu.add(agList(s,j));
            end
        end
    end
    chk = ismember(0,visitedV);
    if (chk == 0)
        break;
    end
end
fprintf('independent components: %d\n',i);
%}

%% find the largest clique

meanVal = mean(mean(agAM));
weightsMap = exp(-agAM/meanVal);
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
[~,endIdx] = min(tmp);

% sort the eigenVectors and consider only the maximum ones
maxMatchCounter = 2000;
currentMatchCounter = 0;
matchesMade = [];
i = 1;
% -- Get threshold
[minVal,minIdx] = min(eigVector);
eigVector(minIdx) = medVal;
[minVal2] = min(eigVector);
thresholdStopVal = minVal2;
eigVector(minIdx) = minVal;



% -- loop starts here
while(currentMatchCounter < maxMatchCounter || i <= agVertexCount)
    [maxVal,maxIdx] = max(eigVector);
    % check if value is less than median of the eigenvector
    if (maxVal <= thresholdStopVal)
        break;
    end
    % get selected Pair
    selPair_idx1 = agVertexList(maxIdx,1);
    selPair_idx2 = agVertexList(maxIdx,2);
    % check for collisions for selected value
    if (isempty(matchesMade))
        tmp1 = [];
        tmp2 = [];
    else
        tmp1  = find(matchesMade(:,1) == selPair_idx1);
        tmp2 = find(matchesMade(:,2) == selPair_idx2);
    end
    if (isempty(tmp2) && isempty(tmp1))
        % pair doesn't collide, add the pair
        tmpMatch = [selPair_idx1, selPair_idx2];
        matchesMade = [matchesMade;tmpMatch];
        currentMatchCounter = currentMatchCounter + 1;
    end
    % make current max value as min value
    eigVector(maxIdx) = minVal;
    % next loop iterator
    i = i+1;
end


denseAssc = zeros(totalVertexCount1,1);
if (~isempty(matchesMade))
    denseAssc(matchesMade(:,1)) = matchesMade(:,2);
end

end