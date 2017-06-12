function plotDenseMap_twinIndex(twinIndex, vertices1, vertices2, landmarks1, landmarks2, offset_x, offset_y, depth_z)

totalVertexCount1 = size(vertices1,2) + size(landmarks1,2);
step_z1           = depth_z / totalVertexCount1;
totalVertexCount2 = size(vertices2,2) + size(landmarks2,2);
step_z2           = depth_z / totalVertexCount2;


hold on;
for i = 1:length(twinIndex)
    p1z = vertices1(twinIndex(i,1)).id*step_z1;
    p2z = vertices2(twinIndex(i,2)).id*step_z2;
    plot3([vertices1(i).x ((vertices2(twinIndex(i)).x)+offset_x)], [vertices1(i).y ((vertices2(twinIndex(i)).y)+offset_y)], [p1z, p2z],  'magenta');
end

end