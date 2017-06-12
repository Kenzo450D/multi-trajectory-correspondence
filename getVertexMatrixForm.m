function [vertexData] = getVertexMatrixForm(vertices)

vCount = size(vertices,2);

vertexData = zeros(4,vCount);
for i = 1:vCount
    vertexData(1,i) = vertices(i).id;
    vertexData(2,i) = vertices(i).x;
    vertexData(3,i) = vertices(i).y;
    vertexData(4,i) = vertices(i).o;
end
end