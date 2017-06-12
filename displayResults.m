function displayResults(denseAssc, landmarkAssc, landmarks1, landmarks2, vertices1, vertices2)


lA1 = landmarkAssc > 0;
lA2tmp = landmarkAssc(landmarkAssc>0);
lA2 = zeros(size(landmarks2,2));
lA2(lA2tmp) = 1;

% [vertices1,vCount1,landmarks1, lCount1, edges1, eCount1, ledges1, leCount1, dim1] =  readLandmarkG2oFile(g2oOpFile1, 'fc.txt', 1);
% [vertices2,vCount2,landmarks2, lCount2, edges2, eCount2, ledges2, leCount2, dim2] =  readLandmarkG2oFile(g2oOpFile2, 'fc.txt', 1);
depth_z = 50;
x_offset = 500;
y_offset = 0;
color1 = jet(size(vertices1,2));
plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, color1,'red');
color2 = getColor2(color1, size(vertices2,2), denseAssc);
% plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, 'blue','red');
plotGraph(vertices2, landmarks2, lA2, x_offset, y_offset, depth_z, color2,'cyan');
plotDenseMap(denseAssc, vertices1, vertices2, landmarks1, landmarks2, x_offset,y_offset, depth_z);

end