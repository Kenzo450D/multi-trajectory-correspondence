function exportLandmarkG2oFile(outFileName, vData, lData, eData, leData)

%% initialization

% initialize counts
vCount  = size(vData, 2);
eCount  = size(eData, 2);
leCount = size(leData, 2);
lCount  = size(lData, 2);
totalVertexCount = vCount + lCount;

vMap = reverseVertexMap(vData, totalVertexCount);

% set up the file print
fileID = fopen(outFileName, 'w');

%% print each of the vertices and edges of each graph

% print the vertices of the graph

lIdx = 1;
for i = 1:totalVertexCount
    if (vMap(i) == 0)
        fprintf(fileID, 'VERTEX_XY %d %d %d\n',(lData(1,lIdx)-1), lData(2,lIdx), lData(3,lIdx));
        lIdx = lIdx + 1;
    else
        fprintf(fileID, 'VERTEX_SE2 %d %d %d %d\n',(vData(1,vMap(i))-1),vData(2,vMap(i)),vData(3,vMap(i)), vData(4,vMap(i)));
    end
end

% printf the edges of the first graph
for i = 1:eCount
    fprintf(fileID, 'EDGE_SE2');
    fprintf(fileID, ' %d', (eData(1,i) - 1));
    fprintf(fileID, ' %d', (eData(2,i) - 1));
    for j = 3:11
        fprintf(fileID, ' %f',eData(j,i));
    end
    fprintf(fileID, '\n');
end

% print the landmark edges of the first graph
for i = 1:leCount
    fprintf(fileID, 'EDGE_SE2_XY');
    fprintf(fileID, ' %d', (leData(1,i) - 1));
    fprintf(fileID, ' %d', (leData(2,i) - 1));
    for j = 3:7
        fprintf(fileID, ' %f',leData(j,i));
    end
    fprintf(fileID, '\n');
end

fclose(fileID);

end