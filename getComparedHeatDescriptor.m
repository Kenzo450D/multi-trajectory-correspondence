function [cmpHd] = getComparedHeatDescriptor(hdv1, hdv2, tScale1, tScale2, poseNeighbourhood1, poseNeighbourhood2 )
%GETCOMPAREDHEADDESCRIPTOR Compares two head descriptors based on the
%information of their neighbourhood sizes, A smaller size means the
%descriptor at the smaller tScale has a higher weightage, while on the
%other hand, the one in a spare locality, would have the higher htScale at
%a higher weightage.
% Input:
%   hdv1: heat descriptor for vertex 1
%   hdv2: heat descriptor for vertex 2
%   tScale1: set of tScales for the first graph
%   tScale2: set of tScales for the second graph
%   poseNeighbourhood1: count of neighbourhood vertices for the particular vertex in graph 1
%   poseNeighbourhood2: count of neighbourhood vertices for the particular vertex in graph 2

% -- Theory:
% Count the number of poses, if the pose count is small, we should give
% more importance to the low t_scale information of the heat kernel
% descriptors for the vertex comparisons. If the pose count is higher, then
% more information for the higher comparisons.

%% increase parameters
increaseParam1 = 10; % multiplier for the tScale chosen
increaseParam2 = 2; % multiplier for the neighbours (1 hop) of the tScale chosen
baseParam = 0.9; % multiplier for the rest of the tScales
ntScales = length(tScale1);

%% checkPoseNeighbourhood with tScale


tmpDiff = abs(tScale1 - poseNeighbourhood1);
[~,tScale1Idx] = min(tmpDiff);

% tmpDiff = abs(tScale2 - poseNeighbourhood2);
% [~,tScale2Idx] = min(tmpDiff);
tScale2Idx = tScale1Idx;

%% make scaling parameter
scaleVector1 = ones(size(tScale1)) * baseParam;
scaleVector1(tScale1Idx) =  increaseParam1;
if (tScale1Idx == ntScales)
    scaleVector1(tScale1Idx - 1) = increaseParam2;
elseif (tScale1Idx == 1)
    scaleVector1(tScale1Idx + 1) = increaseParam2;
else
    scaleVector1(tScale1Idx + 1) = increaseParam2;
    scaleVector1(tScale1Idx - 1) = increaseParam2;
end

scaleVector2 = ones(size(tScale2)) * baseParam;
scaleVector2(tScale2Idx) =  increaseParam1;
if (tScale2Idx == ntScales)
    scaleVector2(tScale2Idx - 1) = increaseParam2;
elseif (tScale2Idx == 1)
    scaleVector2(tScale2Idx + 1) = increaseParam2;
else
    scaleVector2(tScale2Idx + 1) = increaseParam2;
    scaleVector2(tScale2Idx - 1) = increaseParam2;
end

%% scale the heat descriptors
% As we declared vertex descriptors, we have:
% vertexDescriptors = zeros(ntScales, nlandAssc, vCount);
% so a vertex descriptor for a single vertex would be ntScales x nlandAssc

for i = 1:ntScales
    hdv1(i,:) = hdv1(i,:) * scaleVector1(i);
    hdv2(i,:) = hdv2(i,:) * scaleVector2(i);
end

%% make the norm of the difference of head descriptors

cmpHd = norm(hdv1 - hdv2);

end