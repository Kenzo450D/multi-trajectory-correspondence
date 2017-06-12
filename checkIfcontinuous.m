function checkIfcontinuous(edges, vertices)
%CHECKIFCONTINUOUS Checks if the edges are continuous with all vertices

eCount = size(edges,2);
vCount = size(vertices,2);
vCheck = zeros(vCount,1);
for i = 1:eCount
    v2 = edges(i).v2;
    v1 = edges(i).v1;
    if (v2 - v1) == 1
        vCheck(v1) = 1;
    end
end
zeroVal = find(vCheck == 0);
for i = 1:length(zeroVal)
    fprintf(1,'%d\n',zeroVal(i));
end


end