function [distDiff] = getDistanceDiffVertex(vData1, vData2, v1_1Idx, v1_2Idx,v2_1Idx, v2_2Idx)
%GETDISTANCEDIFFVERTEX Returns the difference of the distances in the
%individual graphs.
%Input:
%  vData1: matrix containing information on vertices of graph1
%  vData2: matrix containing information on vertices of graph2
%  v1_1Idx: index of graph1, in vertex 1 of adjacency graph
%  v1_2Idx: index of graph2, in vertex 1 of adjacency graph
%  v2_1Idx: index of graph1, in vertex 2 of adjacency graph
%  v2_2Idx: index of graph2, in vertex 2 of adjacency graph
%Output:
%  distDiff: difference of distance, single value, floating point

% vertex1: graph 1
v1_1x = vData1(2,v1_1Idx);
v1_1y = vData1(3,v1_1Idx);
% vertex1: graph 2
v1_2x = vData2(2,v1_2Idx);
v1_2y = vData2(3,v1_2Idx);
% vertex2: graph 1
v2_1x = vData1(2,v2_1Idx);
v2_1y = vData1(3,v2_1Idx);
% vertex1: graph 2
v2_2x = vData2(2,v2_2Idx);
v2_2y = vData2(3,v2_2Idx);

% distn in graph 1
diff_x = v1_1x - v2_1x;
diff_y = v1_1y - v2_1y;
distnG1 = sqrt(diff_x*diff_x + diff_y*diff_y);
% distn in graph 2
diff_x = v1_2x - v2_2x;
diff_y = v1_2y - v2_2y;
distnG2 = sqrt(diff_x*diff_x + diff_y*diff_y);
% difference
distDiff = distnG1 - distnG2;
end