function [A] = makeIncidenceMatrix(vertices,vCount,landmarks, lCount, edges,eCount, ledges, leCount)
%MAKEADJACENCYMATRIX Makes the adjacency matrix off a landmark graph
%   Makes constanct edge weights based matrix.
%   If odometry edge, then the weight is 0.92, if it's landmark to pose edge,
%   then the weight is 0.7

totaleCount = eCount + leCount;
totalvCount = vCount + lCount;
% -- convert landmarks to an array
lArray = zeros(3,lCount);
for i = 1:lCount
    lArray(1,i) = landmarks(i).id;
    lArray(2,i) = landmarks(i).x;
    lArray(3,i) = landmarks(i).y;
end

edgeWeights = getEdgeWeights(edges, ledges);

%{
~*~*~*If required to modify weights based on distance, add it here*~*~*~
%}

% -- create the incidence matrix
A = zeros(totalvCount, totaleCount);

% -- form the incidence Matrix
for i = 1:eCount
    A(edges(i).v1, i) = -edgeWeights(i);
    A(edges(i).v2, i) = edgeWeights(i);
end
for i = eCount+1:totaleCount
    A(ledges(i - eCount).v1,i) = -edgeWeights(i);
    A(ledges(i - eCount).v2,i) = edgeWeights(i);
end

