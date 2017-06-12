function addSyntheticLandmarkKitti(opFile, outFile)

%% init params
landmarkLookupRange = 20;
landmarkExistRange = 20;


%% Read a kitti g2o file
fileName = opFile;
%fileName = '../input/kitti00op.g2o';

[vertices,vCount,landmarks, lCount, edges,eCount, ledges, leCount, dim] =  readLandmarkG2oFile(fileName , 'fc.txt', 1);

% -- remove loop closures
edges = removeLoopClosures(edges);
eCount = size(edges,2);

% -- convert data to matrix form
[vData, lData, eData, leData] = getMatrixForm(vertices, landmarks, edges, ledges);
eDataOld = eData;

%% find indexes of loop closure edges
lcEdgeIdx = (abs(eData(1,:) - eData(2,:))>1);

%% consider random (X-Y) points as landmarks

% -- calculate average pose size
avgPseStep = getAveragePoseStep(vData);

% -- find the min and max of the vertices
minX = min(vData(2,:));
maxX = max(vData(2,:));
minX = minX(1); % for multiple minimas
maxX = maxX(1); % for multiple maximas
minY = min(vData(3 ,:));
maxY = max(vData(3,:));
minY = minY(1); % for multiple minimas
maxY = maxY(1); % for multiple maximas
meanX = (maxX + minX) / 2;
meanY = (maxY + minY) / 2;
diffX = (maxX - minX);
diffY = (maxY - minY);

% -- randomly generate 100 landmarks
generatedlCount = 1300;
randomPointsX = ((rand(generatedlCount,1) - 0.5)*diffX*2);
randomPointsY = ((rand(generatedlCount,1) - 0.5)*diffY*2);
randomPoints = [randomPointsX, randomPointsY];

% -- remove landmarks which are too close to one another
for i = 1:length(randomPoints)
    for j = 1:length(randomPoints)
        if (i~= j && i < length(randomPoints) && j < length(randomPoints))
            distn = sqrt((randomPoints(i,1) - randomPoints(j,1))^2 + (randomPoints(i,2) - randomPoints(j,2))^2);
            if (distn < 10)
                randomPoints(j,:) = [];
                j = j - 1;
            end
        end
    end
end
generatedlCount = length(randomPoints);

keepLandmark = zeros(generatedlCount,1);
% -- check whether the landmarks are within delta distance from the trajectory, if not, they should be removed.
for i = 1:vCount
    for j = 1:generatedlCount
        lpx = vData(2,i) - randomPoints(j,1);
        lpy = vData(3,i) - randomPoints(j,2);
        lpDist = sqrt(lpx*lpx + lpy*lpy);
        if (lpDist < landmarkExistRange * avgPseStep)
            keepLandmark(j) = keepLandmark(j) + 1;
        end
    end
end

% -- remove the un-needed landmarks
rmLandmarkIdx = find(keepLandmark == 0);
randomPoints(rmLandmarkIdx,:) = [];

%% formulate lData matrix
lCount = length(randomPoints);
lData = zeros(3,lCount);

% -- initialize the landmark pose edge matrix
leData = zeros(7,sum(keepLandmark));
leCount = 0;

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

%% edit the edge list based on the updated vertex information
for i = 1:eCount
    eData(1,i) = vData(1,eData(1,i));
    eData(2,i) = vData(1,eData(2,i));
end

%% edit the ledge list based on the updated vertex information
for i = 1:leCount
    leData(1,i) = vData(1,leData(1,i));
    leData(2,i) = landmarkIdx(leData(2,i));
end


%% plot all of the information in a graph
%figure('name','optimised graph');
%color1 = jet(size(vData,2));
%plotGraph(vData, lData, color1);

%% checks for consistency of vertex and landmark poses
% -- check whether there are no repeating lData index in vData

% -- check for repeat in landmarks
for i = 1:lCount
    for j = i+1:lCount
        if (lData(1,i) == lData(1,j))
            pause(1);
            flagAllOk = 0;
            error('lData has repeating indices');
        end
    end
end


flagAllOk = 1;
for i = 1:lCount
    searchResult = find(vData(1,:) == lData(1,i));
    if (~isempty(searchResult))
        fprintf(1,'ERROR: lData index: %d matches with vData index: %d',i,searchResult(1));
        flagAllOk = 0;
    end
end
if (flagAllOk == 1)
    fprintf(1,'All vertex and landmark indices okay!\n');
end

% -- check whether old eData is retained in the new eData
flagAllOk = 1;
for i = 1:eCount
    v1_old = eDataOld(1,i);
    v2_old = eDataOld(2,i);
    v1_new = vData(1,v1_old);
    v2_new = vData(1,v2_old);
    v1_obs = eData(1,i);
    v2_obs = eData(2,i);
    if ((v1_new == v1_obs) && (v2_new == v2_obs))
        continue;
    else
        flagAllOk = 0;
        error('Error in edge-data data transform\n');
    end
end
if (flagAllOk == 1)
    fprintf(1,'All edges are okay!\n');
end

% -- check whether leData has all of it's v2(2nd pose index) as landmark
% index
flagAllOk = 1;
for i = 1:leCount
    v1 = leData(1,i);
    v2 = leData(2,i);
    vSearchResult = find(vData(1,:) == v1);
    lSearchResult = find(lData(1,:) == v2);
    if(isempty(vSearchResult))
        flagAllOk = 0;
        error('Error in ledges\n');
    end
    if (isempty(lSearchResult))
        flagAllOk = 0;
        error('Error in ledges\n');
    end
end

if (flagAllOk == 1)
    fprintf(1,'All landmark-edges are okay!\n');
end

%% export file to g2o
exportLandmarkG2oFile(outFile, vData, lData, eData, leData);

end