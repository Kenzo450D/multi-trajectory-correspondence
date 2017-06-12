function exportTwinG2oFile(outFileName, vertices1, landmarks1, edges1, ledges1, vertices2, landmarks2, edges2, ledges2, denseMap)

%% initialization
% initialize offset
offset_x = 300;
offset_y = 0;

% initialize counts
vCount1  = size(vertices1, 2);
vCount2  = size(vertices2, 2);
eCount1  = size(edges1, 2);
eCount2  = size(edges2, 2);
leCount1 = size(ledges1, 2);
leCount2 = size(ledges2, 2);
lCount1  = size(landmarks1, 2);
lCount2  = size(landmarks2, 2);
totalVertexCount1 = vCount1 + lCount1;
totalVertexCount2 = vCount2 + lCount2;

% get matrix forms
[vData1, lData1, eData1, leData1] = getMatrixForm(vertices1, landmarks1, edges1, ledges1);
[vData2, lData2, eData2, leData2] = getMatrixForm(vertices2, landmarks2, edges2, ledges2);
vMap1 = reverseVertexMap(vData1, totalVertexCount1);
vMap2 = reverseVertexMap(vData2, totalVertexCount2);

% set up the file print
fileID = fopen(outFileName, 'w');

%% reposition vertices of graph 2 based on the offset
vData2(2,:) = vData2(2,:) + offset_x;
vData2(3,:) = vData2(3,:) + offset_y;
lData2(2,:) = lData2(2,:) + offset_x;
lData2(3,:) = lData2(3,:) + offset_y;

%% print each of the vertices and edges of each graph

% print the vertices of the first graph

lIdx = 1;
for i = 1:totalVertexCount1
    if (vMap1(i) == 0)
        fprintf(fileID, 'VERTEX_XY %d %d %d\n',(lData1(1,lIdx)-1), lData1(2,lIdx), lData1(3,lIdx));
        lIdx = lIdx + 1;
    else
        fprintf(fileID, 'VERTEX_SE2 %d %d %d %d\n',(vData1(1,vMap1(i))-1),vData1(2,vMap1(i)),vData1(3,vMap1(i)), vData1(4,vMap1(i)));
    end
end
lastVIdx1 = (vData1(1,end)-1); % we have increased each index by 1 for matlab

% print the vertices of the second graph

lIdx = 1;
for i = 1:totalVertexCount2
    if (vMap2(i) == 0)
        fprintf(fileID, 'VERTEX_XY %d %d %d\n',(lData2(1,lIdx)+lastVIdx1), lData2(2,lIdx), lData2(3,lIdx));
        lIdx = lIdx + 1;
    else
        fprintf(fileID, 'VERTEX_SE2 %d %d %d %d\n',(vData2(1,vMap2(i))+lastVIdx1),vData2(2,vMap2(i)),vData2(3,vMap2(i)), vData2(4,vMap2(i)));
    end
end

% printf the edges of the first graph
for i = 1:eCount1
    fprintf(fileID, 'EDGE_SE2');
    fprintf(fileID, ' %d', (eData1(1,i) - 1));
    fprintf(fileID, ' %d', (eData1(2,i) - 1));
    for j = 3:11
        fprintf(fileID, ' %d',eData1(j,i));
    end
    fprintf(fileID, '\n');
end

% print the landmark edges of the first graph
for i = 1:leCount1
    fprintf(fileID, 'EDGE_SE2_XY');
    fprintf(fileID, ' %d', (leData1(1,i) - 1));
    fprintf(fileID, ' %d', (leData1(2,i) - 1));
    for j = 3:7
        fprintf(fileID, ' %d',leData1(j,i));
    end
    fprintf(fileID, '\n');
end


% print the edges of the second graph
for i = 1:eCount2
    fprintf(fileID, 'EDGE_SE2');
    fprintf(fileID, ' %d', (eData2(1,i) + lastVIdx1));
    fprintf(fileID, ' %d', (eData2(2,i) + lastVIdx1));
    for j = 3:11
    fprintf(fileID, ' %d',eData2(j,i));
    end
    fprintf(fileID, '\n');
end

% print the landmark edges of the second graph
for i = 1:leCount2
    fprintf(fileID, 'EDGE_SE2_XY');
    fprintf(fileID, ' %d', (leData2(1,i) + lastVIdx1));
    fprintf(fileID, ' %d', (leData2(2,i) + lastVIdx1));
    for j = 3:7
        fprintf(fileID, ' %d',leData2(j,i));
    end
    fprintf(fileID, '\n');
end


%% print the association edges of two graphs

for i = 1:length(denseMap)
    if (denseMap(i) ~= 0)
        % generate and print the vertex indices
        v1Id = i-1;
        v2Id = denseMap(i) + lastVIdx1;
        %fprintf(1, 'EDGE_SE2 %d %d\n',v1Id, v2Id);
        fprintf(fileID, 'EDGE_SE2 %d %d',v1Id, v2Id);
        % generate and print the edge information
        randDisp = rand(1,3) * 0.00001;
        dx = randDisp(1);
        dy = randDisp(2);
        dth = randDisp(3);
        fprintf(fileID, ' %d %d %d',dx, dy, dth);
        % print the information Matrix as well
        for j = 6:11
            fprintf(fileID, ' %d',(eData1(j,1)/100));
        end
        fprintf(fileID, '\n');
    end
end


fclose(fileID);

end