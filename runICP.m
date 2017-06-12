function [Ricp, Ticp, ER, t, newMatchesMade, denseAssc]  = runICP(matchesMade, vertices1, vertices2, landmarks1, landmarks2, vd1, vd2)
%RUNICP Runs ICP on the given initialization on matches. Crops the
%trajectory and works on the required patch.


%% convert to matrix form
vData1 = getVertexMatrixForm(vertices1);
vData2 = getVertexMatrixForm(vertices2);

%% start and end points
tj1v1 = min(matchesMade(:,1));
tj1v2 = max(matchesMade(:,1));
tj2v1 = min(matchesMade(:,2));
tj2v2 = max(matchesMade(:,2));

% -- DEBUG
fprintf(1,'Init : trajectory 1: vertex start: %d\n',tj1v1);
fprintf(1,'Init : trajectory 1: vertex end: %d\n',tj1v2);
fprintf(1,'Init : trajectory 2: vertex start: %d\n',tj2v1);
fprintf(1,'Init : trajectory 2: vertex end: %d\n',tj2v2);
% -- DEBUG END

%% check if the sizes are equal, if not, use "force"
% sz1 = tj1v2 - tj1v1 + 1;
% sz2 = tj2v2 - tj2v1 + 1;
% if (sz1 < sz2)
%     % -- try to increase sz1 by increasing the last point
%     diffn = sz2 - sz1;
%     if ((tj1v2 + diffn) <= size(vData1,2))
%         tj1v2 = tj1v2 + diffn;
%     else
%         canhandleDiff = size(vData1,2) - tj1v2;
%         tj1v2 = size(vData1,2);
%         newdiffn = diffn - canhandleDiff;
%         tj1v1 = tj1v1 - newdiffn;
%     end
% elseif (sz2 < sz1)
%     % -- try to increase sz2 by increasing last point
%     diffn = sz1 - sz2;
%     if ((tj2v2 + diffn) <= size(vData2,2))
%         tj2v2 = tj2v2 + diffn;
%     else
%         canhandleDiff = size(vData2,2) - tj2v2;
%         tj2v2 = size(vData2,2);
%         newdiffn = diffn - canhandleDiff;
%         tj2v1 = tj2v1 - newdiffn;
%     end
% end
% 
% fprintf(1,'Corrected Sizes : trajectory 1: vertex start: %d\n',tj1v1);
% fprintf(1,'Corrected Sizes : trajectory 1: vertex end: %d\n',tj1v2);
% fprintf(1,'Corrected Sizes : trajectory 2: vertex start: %d\n',tj2v1);
% fprintf(1,'Corrected Sizes : trajectory 2: vertex end: %d\n',tj2v2);
    

%% pose information
poseData1 = vData1(2:3,:);
poseData2 = vData2(2:3,:);

%% zero Padding
zp = zeros(1,size(vData1, 2));
poseData1 = [poseData1;zp];
zp = zeros(1,size(vData2, 2));
poseData2 = [poseData2;zp];

%% truncate the poseData to the range declared by start and end points
poseData1 = poseData1(:,tj1v1:tj1v2);
poseData2 = poseData2(:,tj2v1:tj2v2);

%% intialize the first ICP matching from the matches matchesMade
% we cannot loop through matches made, as it would create a problem if the
% first index is changed, else we loop over the index on initMatch
%{
% initMatch = zeros(1,(tj1v2 - tj1v1 + 1));
% for i = 1:size(matchesMade,1)
%     idx = matchesMade(i,1) - tj1v1 + 1;
%     idx2 = matchesMade(i,2) - tj2v1 + 1;
%     initMatch(idx) = idx2;
% end
%}
%
initMatch = zeros(1,(tj1v2 - tj1v1 + 1));
for i = 1:(tj1v2 - tj1v1 + 1)
    idx1 = i + tj1v1 - 1;
    sR = find(matchesMade(:,1) == idx1, 1);
    if (~isempty(sR))
        initMatch(i) = matchesMade(sR, 2) - tj2v1 + 1;
    end
end
%
% we need initMatch to be the opposite mapped as the "match" in the
% algorithm. i.e. We need the indexes of P as mapped as in Q.
% initMatch = zeros(1,(tj2v2 - tj2v1 + 1));
% for i = 1:(tj2v2 - tj2v1 + 1)
%     idx2 = i + tj2v1 - 1;
%     sR = find(matchesMade(:,2) == idx2, 1);
%     if (~isempty(sR))
%         initMatch(i) = matchesMade(sR, 1) - tj1v1 + 1;
%     end
% end


%% run the ICP algorithm
% We run the algorithm only on the associated data points (matchesMade),
% find out the RT necessary, apply that on the entire dataset (poseData),
% and then carry out the rest of the icp. 

[R, T] = icpForMatchesMade(matchesMade, vData1, vData2);

% transform the point cloud for vData2, and run ICP on rest.
poseData2 = R * poseData2 + repmat(T, 1, size(poseData2,2));
[Ricp, Ticp, ER, t, match] = icp(poseData1,poseData2,15);

% An added step should also be done, which showcases the difference between
% the two kinds of initialization.


%[Ricp, Ticp, ER, t, match] = newIcp(poseData1, poseData2, initMatch, 15);

%%  convert match to form of matchesMade
newMatchesMade = zeros(size(match,2),2);
newMatchesMade(:,2) = 1:(tj2v2 - tj2v1 + 1);
newMatchesMade(:,2) = newMatchesMade(:,2) + tj2v1 - 1; % the '1' is to adjust for the matlab indexing (starts from 1 and not 0)
newMatchesMade(:,1) = match + tj1v1 - 1;

%{
% newMatchesMade = zeros((tj1v2 - tj1v1 + 1),2);
% newMatchesMade(:,1) = 1:(tj1v2 - tj1v1 + 1);
% newMatchesMade(:,1) = newMatchesMade(:,1) + tj1v1 - 1; 
% newMatchesMade(:,2) = match + tj2v1 - 1;
%}

%% final check for vertex descriptors in the associations made

mmc = length(newMatchesMade); % matches made count
i = 1;
descr_threshold = 1.6;
while(i<=mmc)
    % -- extract vertex descriptors
    hdv1   = vd1(:,:,newMatchesMade(i,1));
    hdv2   = vd2(:,:,newMatchesMade(i,2));
    % -- compare vertex descriptors
    cmpHd = norm(hdv1 - hdv2);
    % -- compare with threshold
    if (cmpHd > descr_threshold)
        newMatchesMade(i,:) = [];
        i = i-1;
        mmc = length(newMatchesMade);
    end
    % -- loop
    i = i + 1;
end

%% calculate dense matching as is required by plot function
totalVertexCount1 = size(vData1,2) + size(landmarks1,2);
denseAssc = zeros(totalVertexCount1,1);
if (~isempty(newMatchesMade))
    denseAssc(vData1(1,newMatchesMade(:,1))) = vData2(1,newMatchesMade(:,2));
end

end
