function [denseAssc] = getDenseAssociationsAG(distn1, distn2, vd1, vd2, landmarkAssc, vertices1, vertices2, landmarks1, landmarks2, edges1, edges2, ledges1, ledges2, poseAssc)
%GETDENSEASSOCIATIONSAG Get a list of dense associations given a list of
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
neighbourLimit = 5;

% - debug
fprintf(1,'In dense Association function:\n');
fprintf(1,'Graph1: \n');
fprintf(1,'    vCount: %d\n',vCount1);
fprintf(1,'    lCount: %d\n', lCount1);
fprintf(1,'Graph2: \n');
fprintf(1,'    vCount: %d\n',vCount2);
fprintf(1,'    lCount: %d\n', lCount2);

%% Init the correspondence graph
cgVertex = zeros(vCount1,vCount2);
%cg = zeros(vCount1*vCount2, vCount1*vCount2);
cmpHds = zeros(vCount1, vCount2);


%% convert vertex and edge data to matrix
[vData1, lData1, eData1, leData1] = getMatrixForm(vertices1, landmarks1, edges1, ledges1);
[vData2, lData2, eData2, leData2] = getMatrixForm(vertices2, landmarks2, edges2, ledges2);

%% create reverse vertex map
vMap1 = reverseVertexMap(vData1, totalVertexCount1);
vMap2 = reverseVertexMap(vData2, totalVertexCount2);

%% make the correspodence graph
%{
tic;
for i = 1:vCount1
    hdv1 = vd1(:,:,i);
    for j = 1:vCount2
        % compare the two vertices; if similar, add the vertex to the graph
        % cgVertex
        % -- compare heat descriptors
        hdv2 = vd2(:,:,j);
        cmpHd = norm(hdv1 - hdv2);
        cmpHds(i, j) = cmpHd;
        if (cmpHd < 1)
            cgVertex(i,j) = 1;
        end
        % compare the existing vertex, and find if they have a match,
        % if yes, create an edge based on the match, to the cg.

    end
end
tt = toc;
fprintf(1,'Time taken to compare heat descriptors: %f',tt);
pause 1;
%}
%% only compare the neighbours of the computers.
landmarkIdx = (1:length(lData1));
[vn1] = getNeighbours(leData1, lData1,landmarkIdx);
landmarkIdx = (1:length(lData2));
[vn2] = getNeighbours(leData2, lData2, landmarkIdx);

vc1 = length(vn1);
vc2 = length(vn2);

% NOTE: Even comparing 568 and 2588 vertex combinations is too much. We
% consider a uniform sampling of 500 vertices in each map, and compare only
% among them.

%{


cg = sparse(vc1 * vc2, vc1*vc2);
agVertices = [];

for i = 1:vc1
    hdv1 = vd1(:,:,vMap1(vn1(i)));
    for j = 1:vc2
        % compare the two vertices; if similar, add the vertex to the graph
        % cgVertex
        % -- compare heat descriptors
        %fprintf(1,'j: %d,\tvn2(j): %d\t',j,vn2(j));
        %fprintf(1,'vMap2(vn2(j))\n',vMap2(vn2(j)));
        hdv2 = vd2(:,:,vMap2(vn2(j)));
        cmpHd = norm(hdv1 - hdv2);
        cmpHds(vMap1(vn1(i)), vMap2(vn2(j))) = cmpHd;
        if (cmpHd < 1)
            tmpVal = [vMap1(vn1(i)), vMap2(vn2(j))];
            agVertices = [agVertices;tmpVal];
            cgVertex(vMap1(vn1(i)), vMap2(vn2(j))) = 1;
        end
        % compare the existing vertex, and find if they have a match,
        % if yes, create an edge based on the match, to the cg.

    end
end

for i = 1:vc1
    for j = 1:vc2
        if (cgVertex(vMap1(vn1(i)), vMap2(vn2(j))) == 1)
            % -- 
        end
    end
end

%}

%% compare on uniform sampling

% here we consider only the 10th vertex in each of the vertex lists
sampling_rate = 4;

usvCount1 = floor(size(vData1,2)/sampling_rate); %uniformSamplingVertexCount1
usvCount2 = floor(size(vData2,2)/sampling_rate); %uniformSamplingVertexCount2

% ---- make sparse matrix to make their adjacency graph

agVertexList = [];
%cmpHds = ones(usvCount1, usvCount2)*100;
% -- form the vertices
tic;
descr_threshold = 0.05;
for i = 1:usvCount1
    index1 = (i-1)*sampling_rate + 1;
   % v1Idx  = vData1(1,index1);
    hdv1   = vd1(:,:,index1);
    agVertexTmpList = [];
    agIndexTmpList = [];
    cmpHds_arr = [];
    %fprintf(1,'i = %d out of %d\n',i,usvCount1);
    for j = 1:usvCount2
        index2 = (j-1)*sampling_rate + 1;
       % v2Idx = vData2(1,index2);
        hdv2   = vd2(:,:,index2);
        if (norm(hdv1) < 1e-6 || norm(hdv2) < 1e-6)
            continue;
        end
        cmpHd = norm(hdv1 - hdv2); % compare between two heat descriptor for vertices
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

%% only select the top 5 nearest neighbours (in terms of cmpHds) for each vertex

%{
cmpHdsBkp = cmpHds;
agInitList = zeros(agVertexCount,5);

for i = 1:usvCount1
    maxVal = max(cmpHds(i,:));
    for k = 1:5
        [~,minIdx] = min(cmpHds(i,:));
        if (minIdx == i)
            cmpHds(i,minIdx) = maxVal;
            [~,minIdx] = min(cmpHds(i,:));
        end
        agInitList(i,k) = minIdx*10 +1;
        cmpHds(i,minIdx) = maxVal;
    end
end
%}

tic;
agList = zeros(agVertexCount,agVertexCount); %adjacency List
agAM = zeros(agVertexCount,agVertexCount); %adjacency matrix
%agLogical = logical(agVertexCount,
agListEndIndexes = ones(agVertexCount,1);
% get pose step differences in each graph
[eucDiff1] = getPoseStepDiffCombinations(vData1);
[eucDiff2] = getPoseStepDiffCombinations(vData2);
for i = 1:agVertexCount
    v1_1Idx = agVertexList(i,1); % vertex of graph1, in vertex1
    v1_2Idx = agVertexList(i,2); % vertex of graph2, in vertex1
    for j = 1:agVertexCount
        if (i ~= j)
            v2_1Idx = agVertexList(j,1); % vertex of graph1, in vertex2
            v2_2Idx = agVertexList(j,2); % vertex of graph2, in vertex2
            %tic;
            distDiff1 = abs(eucDiff1(v1_1Idx, v2_1Idx) - eucDiff2(v1_2Idx,v2_2Idx));
            %toc;
            %tic;
            %distDiff2 = getDistanceDiffVertex(vData1, vData2, v1_1Idx, v1_2Idx,v2_1Idx, v2_2Idx);
            %toc;
            agAM(i,j) = distDiff1;
            if (distDiff1 < 0.001)
                agList(i,agListEndIndexes(i)) = j;
                agListEndIndexes(i) = agListEndIndexes(i) + 1;
            end
        end 
    end
end
fprintf(1,'Time taken to create edges: %d\n',toc);

%% plot all the correspondences
%{
lA1 = landmarkAssc > 0;
lA2tmp = landmarkAssc(landmarkAssc>0);
lA2 = zeros(size(landmarks2,2));
lA2(lA2tmp) = 1;
depth_z = 50;
x_offset = 500;
y_offset = 0;
plotGraph(vertices1, landmarks1, lA1, edges1, ledges1, 0, 0, depth_z, 'blue','red');
plotGraph(vertices2, landmarks2, lA2, edges2, ledges2, x_offset, y_offset, depth_z, 'green','cyan');
plotDenseMap_twinIndex(agVertexList, vertices1, vertices2, landmarks1, landmarks2, x_offset,y_offset, depth_z);
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
    if (maxVal < thresholdStopVal)
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
denseAssc(matchesMade(:,1)) = matchesMade(:,2);

%% plot the data here..

%{
lA1 = landmarkAssc > 0;
lA2tmp = landmarkAssc(landmarkAssc>0);
lA2 = zeros(size(landmarks2,2));
lA2(lA2tmp) = 1;

% [vertices1,vCount1,landmarks1, lCount1, edges1, eCount1, ledges1, leCount1, dim1] =  readLandmarkG2oFile(g2oOpFile1, 'fc.txt', 1);
% [vertices2,vCount2,landmarks2, lCount2, edges2, eCount2, ledges2, leCount2, dim2] =  readLandmarkG2oFile(g2oOpFile2, 'fc.txt', 1);
depth_z = 50;
x_offset = 500;
y_offset = 0;
plotGraph(vertices1, landmarks1, lA1, edges1, ledges1, 0, 0, depth_z, 'blue','red');
plotGraph(vertices2, landmarks2, lA2, edges2, ledges2, x_offset, y_offset, depth_z, 'green','cyan');
plotDenseMap_twinIndex(denseAssc, vertices1, vertices2, landmarks1, landmarks2, x_offset,y_offset, depth_z);
%}

%{

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

%}
end