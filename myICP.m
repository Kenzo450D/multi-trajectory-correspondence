function myICP(denseAssc, vertices1, vertices2, gtPoseAssc, iterations)
%MYICP A sample code to check icp performance
% -- Add noise
[vertices1] = addNoiseToTrajectoryPoses(vertices1,2);
[vertices2] = addNoiseToTrajectoryPoses(vertices2,2);
% -- convert to matrix
[vData1] = getVertexMatrixForm(vertices1);
[vData2] = getVertexMatrixForm(vertices2);
vMap1 = reverseVertexMap(vData1);
vMap2 = reverseVertexMap(vData2);

v1Poses = find(denseAssc>0);
v2Poses = denseAssc(denseAssc>0);
v1MappedPoses = vMap1(v1Poses);
v2MappedPoses = vMap2(v2Poses);

% -- ground Truth ICP
v1gt = find(gtPoseAssc>0);
v2gt = gtPoseAssc(gtPoseAssc>0);

% -- send to ICP
p = vData1(2:3,v1gt);
q = vData2(2:3,v2gt);

q = q(:,1:10);
p = p(:,1:10);
addZero = zeros(1,size(q,2));
q = [q;addZero];
p = [p;addZero];
initMatch = [1,0,3,0,5,0,7,0,9,0];
%[Ricp Ticp ER t match] = newIcp(p, q, 15);
[Ricp, Ticp, ER, t, match] = newIcp(p, q, [], 15);
fprintf(1,'Work without initial matching');
disp(ER);
fprintf(1,'Matching:');
disp(match);
fprintf(1,'Work with initial matching');
[Ricp, Ticp, ER, t, match] = newIcp(p, q, initMatch, 15);
disp(ER);
fprintf(1,'Matching(with init mapping):');
disp(match);
end