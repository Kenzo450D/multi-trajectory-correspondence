function plotGraphWithoutLandmarks(vertices, offset_x, offset_y, depth_z, color1)

totalVertexCount = size(vertices,2);
step_z           = depth_z / totalVertexCount;

hold on;
vData = getVertexMatrixForm(vertices);
vData(2,:) = vData(2,:) + offset_x;
vData(3,:) = vData(3,:) + offset_y;
offset_z = linspace(0,depth_z,totalVertexCount);
vzSteps = offset_z;
scatter3(vData(2,:), vData(3,:), vzSteps,5,color1,'filled');
hold off;
end