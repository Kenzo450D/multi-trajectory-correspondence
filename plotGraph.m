function plotGraph(vertices, landmarks, landmarkAssc, offset_x, offset_y, depth_z, color1, color2)

totalVertexCount = size(vertices,2) + size(landmarks,2);
step_z           = depth_z / totalVertexCount;

hold on;
vData = getVertexMatrixForm(vertices);
vData(2,:) = vData(2,:) + offset_x;
vData(3,:) = vData(3,:) + offset_y;
offset_z = linspace(0,depth_z,totalVertexCount);
vzSteps = offset_z(vData(1,:));
scatter3(vData(2,:), vData(3,:), vzSteps,5,color1,'filled');

for i = 1:size(landmarks,2)
    if(landmarkAssc(i) == 1)
        scatter3((landmarks(i).x+offset_x), (landmarks(i).y+offset_y), ((landmarks(i).id)*step_z), 20,'black','filled');
    else
        scatter3((landmarks(i).x+offset_x), (landmarks(i).y+offset_y), ((landmarks(i).id)*step_z), 10,color2,'filled');
    end
end
%{
for i= 1:size(edges,2)
    plot([(vertices(edges(i).v1).x+offset_x) (vertices(edges(i).v2).x+offset_x)],[(vertices(edges(i).v1).y+offset_y) (vertices(edges(i).v2).y+offset_y)]);
end
for i = 1:size(ledges,2)
    plot([(vertices(edges(i).v1).x+offset_x) (vertices(edges(i).v2).x+offset_x)],[(landmarks(edges(i).v1).y+offset_y) (landmarks(edges(i).v2).y+offset_y)]);
end
%}
hold off;
end