function [matchesMade, denseAssc] = getDenseAssociationsLocal(distn1, distn2, vd1, vd2, landmarkAssc, vertices1, vertices2, landmarks1, landmarks2, edges1, edges2, ledges1, ledges2, poseAssc, tScale1, tScale2)
%GETDENSEASSOCIATIONSLOCAL Get a list of dense associations given a list of
%landmark Associations, and the vertex and edge informations. This uses
%Associate Graph, and also the fact that a local matching is only made, with
%neighbours of landmark to one another
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
landmarkNLimit = 50; %neighbours around landmarks which are considered
descr_threshold = 0.05;

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

%% only compare the neighbours of the landmarks
landmarkIdx = (1:length(lData1));
[vn1] = getNeighbours(leData1, lData1,landmarkIdx);
landmarkIdx = (1:length(lData2));
[vn2] = getNeighbours(leData2, lData2, landmarkIdx);

vc1 = length(vn1);
vc2 = length(vn2);

%% compare with local subsets of neighbours
% We have to take a local subset around every landmark. Now we know that there
% exists landmarks close to one another, and we can't always make the right
% choices in selecting only the 50 neighbours around one landmark as there
% would be overlaps in selections.

% Step 1:  Choose landmarkNLimit around each landmark, store it in structures
%      1.a It would be lesser when a trajectory ends at one side
% Step 2:  Merge the ones which overlap with the other ones.
% Step 3:  Calculate associate Graph function, to get dense asssociation for
%          two subgraphs.

l1Assc = landmarkAssc>0;
tmpVar = [1:size(landmarks1,2)];
l1Assc = tmpVar(l1Assc);
[indexEnds1] = getDynamicSubsetNeighbours(vMap1, vData1, lData1, l1Assc, landmarkNLimit);

l2Assc = landmarkAssc>0;
l2Assc = landmarkAssc(l2Assc);
[indexEnds2] = getDynamicSubsetNeighbours(vMap2, vData2, lData2, l2Assc, landmarkNLimit);

% get pose step differences in each graph
[eucDiff1] = getPoseStepDiffCombinations(vData1);
[eucDiff2] = getPoseStepDiffCombinations(vData2);

%% get neighbourhood information

poseNeighbourhoodMap1 = zeros(vCount1, 1);
poseNeighbourhoodMap2 = zeros(vCount2, 1);
for i = 1:(length(l1Assc))
    v1start = indexEnds1(i,1);
    v1end   = indexEnds1(i,2);
    poseDiffn = v1end - v1start;
    for j = v1start:v1end
        poseNeighbourhoodMap1(j) = poseDiffn;
    end
    v2start = indexEnds2(i,1);
    v2end   = indexEnds2(i,2);
    poseDiffn = v2end - v2start;
    for j = v2start:v2end
        poseNeighbourhoodMap2(j) = poseDiffn;
    end
end

%% amplify heat descriptors
[vd1] = amplifyHeatDescriptor(poseNeighbourhoodMap1, vCount1, vd1);
[vd2] = amplifyHeatDescriptor(poseNeighbourhoodMap2, vCount2, vd2);

% -- START DEBUG
%{
% for i = 1:(length(l1Assc))
%     v1start = indexEnds1(i,1);
%     v1end   = indexEnds1(i,2);
%     v2start = indexEnds2(i,1);
%     v2end   = indexEnds2(i,2);
%     diffn = v1end - v1start;
%     fprintf(1,'section %d:\t%d\n',i,diffn);
% end
%}
% -- END DEBUG

%% get the vertex assocations
agAssociations = [];
denseAssc = zeros(totalVertexCount1,1);
for i = 1:(length(l1Assc)-3)
    v1start = indexEnds1(i,1);
    v1end   = indexEnds1(i+3,2);
    v2start = indexEnds2(i,1);
    v2end   = indexEnds2(i+3,2);
    %{
    if ((i < length(l1Assc)) && (i > 1))
        v2start = indexEnds2(i-1,1);
        v2end   = indexEnds2(i+1,2);
    elseif (i == 1)
        v2start = indexEnds2(i,1);
        v2end = indexEnds2(i+1,2);
    elseif (i == length(l1Assc))
        v2start = indexEnds2(i-1,1);
        v2end = indexEnds2(i,2);
    end
    %}
    % -- START DEBUG
    %fprintf(1,'i = %d out of %d\n',i,length(l1Assc));
    %fprintf(1,'........v1Start: %d\tv1end: %d\n........v2start: %d\tv2end: %d\n',v1start,v1end, v2start, v2end);
    %tmpDenseAssc = getDenseAssociationsSubset(vd1, vd2, v1start, v1end, v2start, v2end, eucDiff1, eucDiff2, neighbourLimit, totalVertexCount1, tScale1, tScale2, poseNeighbourhoodMap1, poseNeighbourhoodMap2);
    %denseAssc(zeroIdx) = tmpAssociations(zeroIdx);
    [tmpAssociations] = getAssociationsVertexSubset(vd1, vd2, v1start, v1end, v2start, v2end, neighbourLimit, totalVertexCount1, tScale1, tScale2, poseNeighbourhoodMap1, poseNeighbourhoodMap2);
    % insert the data into denseAssc for the ones which are zero 
    zeroIdx = denseAssc == 0;
    %disp(tmpAssociations);
    %fprintf(1,'adding associations to the total set, next loop!\n');
    agAssociations = [agAssociations;tmpAssociations];
end

%% use association graph to prune outlier associations, and form the final output

[matchesMade] = associationGraphDenseAssocation(agAssociations, eucDiff1, eucDiff2, ...
                                 totalVertexCount1, vertices1, vertices2, ...
                                 vd1, vd2,...
                                 landmarks1, landmarks2, landmarkAssc);
            %associationGraphDenseAssocation(agVertexList, eucDiff1, eucDiff2, totalVertexCount1, vertices1, vertices2, landmarks1, landmarks2, landmarkAssc)

fprintf(1,'Dense Associations made\n');

%% final check for vertex descriptors in the associations made

descr_threshold = 2;
mmc = length(matchesMade); % matches made count
for i = 1:mmc
    hdv1   = vd1(:,:,matchesMade(i,1));
    hdv2   = vd2(:,:,matchesMade(i,2));
    cmpHd = norm(hdv1 - hdv2);
    if (cmpHd > descr_threshold)
        matchesMade(i,:) = [];
        i = i-1;
        mmc = mmc - 1;
    end
end

denseAssc = zeros(totalVertexCount1,1);
if (~isempty(matchesMade))
    denseAssc(vData1(1,matchesMade(:,1))) = vData2(1,matchesMade(:,2));
end


end