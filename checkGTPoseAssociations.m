function checkGTPoseAssociations(l1Assc, l2Assc, vd1, vd2, landmarkAssc, vertices1, vertices2,landmarks1, landmarks2, edges1, edges2, ledges1, ledges2, poseAssc);
%CHECKGTPOSEASSOCIATIONS check how much the vertex descriptors match up to
%each other in case of ground truth. Ideally values should be zero.

%% initializations
vCount1 = size(vertices1, 2);
vCount2 = size(vertices2, 2);
lCount1 = size(landmarks1, 2);
lCount2 = size(landmarks2, 2);
totalVertexCount1 = vCount1 + lCount1;
totalVertexCount2 = vCount2 + lCount2;
neighbourLimit = 5;
landmarkNLimit = 50; %neighbours around landmarks which are considered

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

l1Assc = landmarkAssc>0;
tmpVar = [1:size(landmarks1,2)];
l1Assc = tmpVar(l1Assc);
[indexEnds1] = getDynamicSubsetNeighbours(vMap1, vData1, lData1, l1Assc, landmarkNLimit);

l2Assc = landmarkAssc>0;
l2Assc = landmarkAssc(l2Assc);
[indexEnds2] = getDynamicSubsetNeighbours(vMap2, vData2, lData2, l2Assc, landmarkNLimit);

%% get pose step differences in each graph
[eucDiff1] = getPoseStepDiffCombinations(vData1);
[eucDiff2] = getPoseStepDiffCombinations(vData2);

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

%% compare heat descriptors

cmpHds = 0;
totalCmp = 0;
allCmpHds = zeros(1,size(poseAssc,1));
for i = 1:length(poseAssc)
    if (poseAssc(i) ~= 0)
        fprintf(1,'i = %d\tposeAssc: %d\n',i,poseAssc(i));
        cmpHd = norm(vd1(:,:,i) - vd2(:,:,poseAssc(i)));
        fprintf(1,':::::::: %d\n',cmpHd);
        totalCmp = totalCmp + 1;
        cmpHds = cmpHds + cmpHd;
        allCmpHds(i) = cmpHd;
    end
end
avgCmpHds = cmpHds / totalCmp;
fprintf(1,'Average compare: %d\n',avgCmpHds);
maxCmpHds = max(allCmpHds);
fprintf(1,'Max Compare: %d\n',maxCmpHds);
minCmpHds = min(allCmpHds);
fprintf(1,'Min Compare: %d\n',minCmpHds);


end