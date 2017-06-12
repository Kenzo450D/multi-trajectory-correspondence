function [edgeWeights] = getEdgeWeights(edges, ledges)

eCount = size(edges,2);
leCount = size(ledges,2);
edgeWeights = zeros(1,eCount + leCount);

edgeWeights(1:eCount) = 0.92;
edgeWeights(eCount+1:eCount+leCount) = 0.8;

end
