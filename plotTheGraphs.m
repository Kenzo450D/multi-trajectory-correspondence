function plotTheGraphs(fileName1, fileName2, landmarkAssc, denseAssc)
%PLOTTHEGRAPHS Plots the graph comparing the two trajectories on a single plot
%  Input:
%   fileName1: .mat fileName for file1 on which to plot
%   fileName2: .mat fileName for file2 on which to plot

% ---- file 1
load(fileName1);
vertices1 = vertices;
landmarks1 = landmarks;
edges1 = edges;
ledges1 = ledges;
% ---- file 2
load(fileName2);
vertices2 = vertices;
landmarks2 = landmarks;
edges2 = edges;
ledges2 = ledges;

% --  visualize data
% ---- create a figure
figure('name','matches found');

% ---- map the landmark matches to each trajectory
lA1 = landmarkAssc > 0;
lA2tmp = landmarkAssc(landmarkAssc>0);
lA2 = zeros(size(landmarks2,2));
lA2(lA2tmp) = 1;

% ---- set the parameters for plot
depth_z = 50;
x_offset = 500;
y_offset = 0;
color1 = jet(size(vertices1,2));
color2 = getColor2(color1, size(vertices2,2), denseAssc);

% ---- plot graph1
plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, color1,'red');

% ---- plot graph2
% plotGraph(vertices1, landmarks1, lA1, 0, 0, depth_z, 'blue','red');
plotGraph(vertices2, landmarks2, lA2, x_offset, y_offset, depth_z, color2,'cyan');

% ---- plot the dense matching between the graphs
plotDenseMap(denseAssc, vertices1, vertices2, landmarks1, landmarks2, x_offset,y_offset, depth_z);

end