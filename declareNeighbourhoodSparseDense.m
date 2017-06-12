function [poseNeighbourhoodType] = declareNeighbourhoodSparseDense(poseNeighbourhoodMap, vCount)
%DECLARENEIGHBOURHOODSPARSEDENSE Declare whether a neighbourhood is sparse
%or dense, based on the mean shift.

% -- check the locality
% 6% - sparse

threshold = 0.06 * vCount;

% we are just concerned about which are ones that are spare, the dense ones
% would just use the direct comparison without amplification of heat
% descriptors.

% so we make poseNeighbourhoodType to be a binary vector, 1 declares
% sparse, and 0 declares otherwise.

poseNeighbourhoodType = zeros(size(poseNeighbourhoodMap));
idx = poseNeighbourhoodMap > threshold;
poseNeighbourhoodType(idx) = 1;

end