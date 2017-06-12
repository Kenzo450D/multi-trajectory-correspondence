function [Ricp, Ticp, ER, t, newMatchesMade]  = runICPOnly(vertices1, vertices2)
%RUNICP Runs ICP on the given initialization on matches. Crops the
%trajectory and works on the required patch.


% -- convert to matrix form
vData1 = getVertexMatrixForm(vertices1);
vData2 = getVertexMatrixForm(vertices2);

% -- start and end points
tj1v1 = 1;
tj1v2 = size(vData1,2);
tj2v1 = 1;
tj2v2 = size(vData2,2);

% -- DEBUG
% fprintf(1,'Init : trajectory 1: vertex start: %d\n',tj1v1);
% fprintf(1,'Init : trajectory 1: vertex end: %d\n',tj1v2);
% fprintf(1,'Init : trajectory 2: vertex start: %d\n',tj2v1);
% fprintf(1,'Init : trajectory 2: vertex end: %d\n',tj2v2);
% -- DEBUG END

% -- check if the sizes are equal, if not, use "force"
% sz1 = tj1v2 - tj1v1 + 1;
% sz2 = tj2v2 - tj2v1 + 1;
% if (sz1 < sz2)
%     % -- reduce the one with lower vertices
%     tj2v2 = tj1v2;
% elseif (sz2 < sz1)
%     tj1v2 = tj2v2;
% end
    
% fprintf(1,'Corrected Sizes : trajectory 1: vertex start: %d\n',tj1v1);
% fprintf(1,'Corrected Sizes : trajectory 1: vertex end: %d\n',tj1v2);
% fprintf(1,'Corrected Sizes : trajectory 2: vertex start: %d\n',tj2v1);
% fprintf(1,'Corrected Sizes : trajectory 2: vertex end: %d\n',tj2v2);
    

% -- pose information
poseData1 = vData1(2:3,:);
poseData2 = vData2(2:3,:);

% -- zero Padding
zp = zeros(1,size(vData1, 2));
poseData1 = [poseData1;zp];
zp = zeros(1,size(vData2, 2));
poseData2 = [poseData2;zp];

% -- truncate the poseData to the range declared by start and end points
poseData1 = poseData1(:,tj1v1:tj1v2);
poseData2 = poseData2(:,tj2v1:tj2v2);

% -- intialize the first ICP matching from the matches matchesMade
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

% -- run the ICP algorithm
[Ricp, Ticp, ER, t, match] = icp(poseData1, poseData2, 15);

% -- convert match to form of matchesMade
newMatchesMade = zeros(size(match,2),2);
newMatchesMade(:,2) = 1:(tj2v2 - tj2v1 + 1);
newMatchesMade(:,2) = newMatchesMade(:,2) + tj2v1 - 1; % the '1' is to adjust for the matlab indexing (starts from 1 and not 0)
newMatchesMade(:,1) = match + tj1v1 - 1;

end
