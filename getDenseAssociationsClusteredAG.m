function [denseAssc] = getDenseAssociationsClusteredAG(distn1, distn2, vd1, vd2, landmarkAssc, vertices1, vertices2, landmarks1, landmarks2, edges1, edges2, ledges1, ledges2, poseAssc)
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
%cgVertex = zeros(vCount1,vCount2);
%cg = zeros(vCount1*vCount2, vCount1*vCount2);
%cmpHds = zeros(vCount1, vCount2);

% The rest of this section is basically code optimisation to make the make
% matlab code run faster, as it's faster to extract data from a matrix than
% a structure.

% --  convert vertex and edge data to matrix
[vData1, lData1, eData1, leData1] = getMatrixForm(vertices1, landmarks1, edges1, ledges1);
[vData2, lData2, eData2, leData2] = getMatrixForm(vertices2, landmarks2, edges2, ledges2);

% --  create reverse vertex map
vMap1 = reverseVertexMap(vData1, totalVertexCount1);
vMap2 = reverseVertexMap(vData2, totalVertexCount2);

%% cluster the graph into sections of high and low density

% Make a variable for every vertex, which shows the sum of the reciprocal
% of the distance from that vertex to all of the landmarks.

vlDist1 = sumOfReciprocatedDistance(vData1,lData1, eData1, leData1);

for i = 1:length(vData1)
    
end

%% compare on uniform sampling

% here we consider only the 10th vertex in each of the vertex lists
sampling_rate = 2;

usvCount1 = floor(size(vData1,2)/sampling_rate); %uniformSamplingVertexCount1
usvCount2 = floor(size(vData2,2)/sampling_rate); %uniformSamplingVertexCount2

% ---- make sparse matrix to make their adjacency graph

agVertexList = [];
%cmpHds = ones(usvCount1, usvCount2)*100;
% -- form the vertices
tic;
descr_threshold = 0.1;
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
        if (norm(hdv1) < 1e-10 || norm(hdv2) < 1e-10)
            continue;
        end
        cmpHd = norm(hdv1 - hdv2);
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
     else         
            disp(sprintf('No Match Found for index: %d in left trajectory',index1));
     end
    %fprintf(1,'length of agVertexList: %d\n',length(agVertexList));
end

agVertexCount = length(agVertexList);
fprintf(1,'Time taken to create vertices: %d\n',toc);

%% only select the top 5 nearest neighbours (in terms of cmpHds) for each vertex


tic;
agList = zeros(agVertexCount,agVertexCount); %adjacency List
agAM = zeros(agVertexCount,agVertexCount); %adjacency matrix
%agLogical = logical(agVertexCount,
agListEndIndexes = ones(agVertexCount,1);
% get euclidean differences in each graph
[eucDiff1] = getPoseDiffCombinations(vData1);
[eucDiff2] = getPoseDiffCombinations(vData2);
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

end