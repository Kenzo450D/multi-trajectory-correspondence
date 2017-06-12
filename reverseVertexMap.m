function [vMap] = reverseVertexMap(vData, totalVertexCount)
%REVERSEVERTEXMAP Reverse Vertex Map
if (exist('totalVertexCount','var'))
    vMap = zeros(1,totalVertexCount);
else
    maxVal = max(vData(1,:));
    vMap = zeros(1,maxVal);
end
for i = 1:length(vData)
    vMap(vData(1,i)) = i;
end
end