function plotDenseMap(denseMap, vertices1, vertices2, landmarks1, landmarks2, offset_x, offset_y, depth_z, skipPlot)

totalVertexCount1 = size(vertices1,2) + size(landmarks1,2);
step_z1           = depth_z / totalVertexCount1;
totalVertexCount2 = size(vertices2,2) + size(landmarks2,2);
step_z2           = depth_z / totalVertexCount2;


%% convert vertex and edge data to matrix
[vData1] = getVertexMatrixForm(vertices1);
[vData2] = getVertexMatrixForm(vertices2);

%% create reverse vertex map
vMap1 = reverseVertexMap(vData1);
vMap2 = reverseVertexMap(vData2);

%% plot the graph
set(gcf,'Color',[1,1,1])
hold on;
for i = 1:length(denseMap)
    if (denseMap(i) ~= 0)
        if(exist('skipPlot','var') && skipPlot ~= 0)
            if (mod(i,skipPlot) ~= 0) 
                continue;
            end
        end
        p1z = vertices1(i).id*step_z1;
        p2z = vertices2(vMap2(denseMap(i))).id*step_z2;
        plot3([vertices1(vMap1(i)).x ((vertices2(vMap2(denseMap(i))).x)+offset_x)], [vertices1(vMap1(i)).y ((vertices2(vMap2(denseMap(i))).y)+offset_y)], [p1z, p2z], 'Color', [0.8,0.8,0.8] );
    end
end

end