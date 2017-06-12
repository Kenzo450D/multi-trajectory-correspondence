function plotDenseMapWithoutLandmarks(matchesMade, vertices1, vertices2, cropVertices1, cropVertices2, offset_x, offset_y, depth_z)

totalVertexCount1 = size(cropVertices1,2);
step_z1           = depth_z / totalVertexCount1;
totalVertexCount2 = size(cropVertices2,2);
step_z2           = depth_z / totalVertexCount2;


%% convert vertex and edge data to matrix
[vData1] = getVertexMatrixForm(vertices1);
[vData2] = getVertexMatrixForm(vertices2);

%% create reverse vertex map
vMap1 = reverseVertexMap(vData1);
vMap2 = reverseVertexMap(vData2);

%% plot the graph
set(gcf,'Color',[1,1,1]);

% vData1_x = vData1(2,matchesMade(:,1));
% vData1_y = vData1(3,matchesMade(:,1));
% 
% vData2_x = vData2(2,matchesMade(:,2)) + offset_x;
% vData2_y = vData2(3,matchesMade(:,2)) + offset_y;
% 
% offset_z1 = linspace(0,depth_z,totalVertexCount1);
% offset_z2 = linspace(0,depth_z,totalVertexCount2);
% 
% vData1_z = offset_z1(matchesMade(:,1) - min(matchesMade(:,1)) + 1);
% vData2_z = offset_z2(matchesMade(:,2) - min(matchesMade(:,2)) + 1);
% 
% plot3([vData1_x vData2_x],[vData1_y vData2_y], [vData1_z vData2_z],'Color', [0.8,0.8,0.8] );

cropPointT1v1 = min(matchesMade(:,1));
cropPointT1v2 = max(matchesMade(:,1));
cropPointT2v1 = min(matchesMade(:,2));
cropPointT2v2 = max(matchesMade(:,2));

hold on;
for i = 1:length(matchesMade)
    p1z = (matchesMade(i,1) - cropPointT1v1 +1)*step_z1;
    p2z = (matchesMade(i,2) - cropPointT2v1 +1)*step_z2;
    plot3([vData1(2,matchesMade(i,1)) ((vData2(2,matchesMade(i,2)))+offset_x)], [ vData1(3,matchesMade(i,1)) (vData2(3,matchesMade(i,2)) + offset_y)], [p1z, p2z], 'Color', [0.8,0.8,0.8] );
%     if(mod(i,100) == 0)
%         pause();
%     end
    %plot3([vertices1(vMap1(i)).x ((vertices2(vMap2(matchesMa(i))).x)+offset_x)], [vertices1(vMap1(i)).y ((vertices2(vMap2(matchesMa(i))).y)+offset_y)], [p1z, p2z], 'Color', [0.8,0.8,0.8] );
end

end