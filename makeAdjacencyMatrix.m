function [A] = makeAdjacencyMatrix(vertices,landmarks, edges, ledges)
%MAKEADJACENCYMATRIX Makes the adjacency matrix off a landmark graph
%   Makes constanct edge weights based matrix.
%   If odometry edge, then the weight is 0.92, if it's landmark to pose edge,
%   then the weight is 0.7

vCount = size(vertices,2);
lCount = size(landmarks,2);
eCount = size(edges,2);
leCount = size(ledges,2);

odometryWeight = 0.92;
loopClosureWeight = 0.7;
landmarkPoseWeight = 0.8;

totalvCount = vCount + lCount;

A = zeros(totalvCount, totalvCount);

for i = 1:eCount
    v1 = edges(i).v1;
    v2 = edges(i).v2;
    if (abs(v1 - v2) == 1 || abs(v1 - v2) == 2) % odometry edge
        A(v1,v2) = odometryWeight;
        A(v2,v1) = odometryWeight;
    else
        A(v1,v2) = loopClosureWeight;
        A(v2,v1) = loopClosureWeight;
    end
end

for i = 1:leCount
    v1 = ledges(i).v1;
    v2 = ledges(i).v2;
    A(v1,v2) = landmarkPoseWeight;
    A(v2,v1) = landmarkPoseWeight;
end
end
