function plotTheCroppedGraphs(opFileName1, opFileName2, newMatchesMade, cropPointT1v1, cropPointT1v2, cropPointT2v1, cropPointT2v2, nameFig)
%PLOTTHEGRAPHS Plots the graph comparing the two trajectories on a single plot
%  Input:
%   opFileName1: .mat fileName for file1 on which to plot
%   opFileName2: .mat fileName for file2 on which to plot

% ---- file 1
load(opFileName1);
vertices1 = vertices;
landmarks1 = landmarks;
edges1 = edges;
ledges1 = ledges;
% ---- file 2
load(opFileName2);
vertices2 = vertices;
landmarks2 = landmarks;
edges2 = edges;
ledges2 = ledges;

% -- crop the data
cropVertices1 = vertices1(cropPointT1v1:cropPointT1v2);
cropVertices2 = vertices2(cropPointT2v1:cropPointT2v2);

% --  visualize data
% ---- create a figure
figure('name',nameFig);

% ---- set the parameters for plot
depth_z = 50;
offset_x = 500;
offset_y = 0;
% color1 = jet(size(vertices1,2));
% color2 = getColor2(color1, size(vertices2,2), denseAssc);
color1 = 'blue';
color2 = 'blue';
% ---- plot graph1
plotGraphWithoutLandmarks(cropVertices1, 0, 0, depth_z, color1);
%plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, color1,'red');

% ---- plot graph2
% plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, 'blue','red');
hold on;
plotGraphWithoutLandmarks(cropVertices2, offset_x, offset_y, depth_z, color2);
%plotGraph(vertices2, landmarks2, lA2, x_offset, y_offset, depth_z, color2,'cyan');

% -- convert matchesMade with respect to the to our cropped points
%newMatchesMade(:,1) = newMatchesMade(:,1) - cropPointT1v1 + 1;
%newMatchesMade(:,2) = newMatchesMade(:,2) - cropPointT2v1 + 1;
hold on;
% ---- plot the dense matching between the graphs
plotDenseMapWithoutLandmarks(newMatchesMade, vertices1, vertices2, cropVertices1, cropVertices2, offset_x, offset_y, depth_z);
hold off;
end