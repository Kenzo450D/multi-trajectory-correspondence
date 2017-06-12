function exportLandmarkG2oGraph(vertices,landmarks, edges, ledges, outputFileName)
% -- export the data to a g2o file

% -- get the size of each of the inputs
vCount = size(vertices,2);
lCount = size(landmarks,2);
eCount = size(edges,2);
leCount = size(ledges,2);


% -- convert the structure data to an array
vertexData = zeros(4,vCount);
for i = 1:vCount
    vertexData(1,i) = vertices(i).id;
    vertexData(2,i) = vertices(i).x;
    vertexData(3,i) = vertices(i).y;
    vertexData(4,i) = vertices(i).o;
end
landmarkData = zeros(3,lCount);
for i = 1:lCount
    landmarkData(1,i) = landmarks(i).id;
    landmarkData(2,i) = landmarks(i).x;
    landmarkData(3,i) = landmarks(i).y;
end
edgeData = zeros(11,eCount);
for i = 1:eCount
    edgeData(1,i) = edges(i).v1;
    edgeData(2,i) = edges(i).v2;
    edgeData(3,i) = edges(i).dx;
    edgeData(4,i) = edges(i).dy;
    edgeData(5,i) = edges(i).dth;
    %fill Covariance Matrix
    edgeData(6,i)  = edges(i).covMatrix(1,1);
    edgeData(7,i)  = edges(i).covMatrix(1,2);
    edgeData(8,i)  = edges(i).covMatrix(1,3);
    edgeData(9,i)  = edges(i).covMatrix(2,2);
    edgeData(10,i) = edges(i).covMatrix(3,2);
    edgeData(11,i) = edges(i).covMatrix(3,3);
end
landmarkEdgeData = zeros(7,leCount);
for i = 1:1:leCount
    landmarkEdgeData(1,i) = ledges(i).v1;
    landmarkEdgeData(2,i) = ledges(i).v2;
    landmarkEdgeData(3,i) = ledges(i).dx;
    landmarkEdgeData(4,i) = ledges(i).dy;
    landmarkEdgeData(5,i) = ledges(i).covMatrix(1,1);
    landmarkEdgeData(6,i) = ledges(i).covMatrix(1,2);
    landmarkEdgeData(7,i) = ledges(i).covMatrix(2,2);
end
    
% -- find the minimum of the vertices, truncate the values accordingly
minVertexId = min(vertexData(1,:));
if (minVertexId > min(landmarkData(1,:)))
    minVertexId = min(landmarkData(1,:));
end

% -- adjust the rest of the edges and vertices accordingly so that the
% minimum is zero

vertexData(1,:) = vertexData(1,:) - minVertexId;
landmarkData(1,:) = landmarkData(1,:) - minVertexId;
edgeData(1,:) = edgeData(1,:) - minVertexId;
edgeData(2,:) = edgeData(2,:) - minVertexId;
landmarkEdgeData(1,:) = landmarkEdgeData(1,:) - minVertexId;
landmarkEdgeData(2,:) = landmarkEdgeData(2,:) - minVertexId;

% -- open the file
outFile = fopen(outputFileName,'w');

% -- print the vertices in the file
lidx = 1;
vidx = 1;
justPrintFlag = 0;
while (true)
    if (lidx <= lCount && vidx <= vCount)
        lidData = landmarkData(1,lidx);
        vidData = vertexData(1,vidx);
    else
        justPrintFlag = 1;
    end
    if ( lidx <= lCount && ( lidData < vidData || justPrintFlag == 1))
        % -- print the landmark vertex
        fprintf(outFile, 'VERTEX_XY %d %f %f\n', landmarkData(1,lidx), landmarkData(2,lidx), landmarkData(3,lidx));
        % -- increment the landmark index to read the next one
        lidx = lidx + 1;
    elseif ( vidx <= vCount && ( lidData > vidData || justPrintFlag == 1))
        % -- print the odometry pose vertex
        fprintf(outFile, 'VERTEX_SE2 %d %f %f %f\n', vertexData(1,vidx), vertexData(2, vidx), vertexData(3,vidx), vertexData(4,vidx));
        % -- increment the vertex index to read the next one
        vidx = vidx+ 1;
    end
    if (lidx > lCount && vidx > vCount)
        break;
    end
end

% -- print the edges in the file
lEIdx = 1;
eIdx = 1;
justPrintFlag = 0;
while(true)
    if (lEIdx <= leCount && eIdx <= eCount)
        lidData = landmarkEdgeData(1,lEIdx);
        eidData = edgeData(1,eIdx);
    else
        justPrintFlag = 1;
    end
    if ( lEIdx <= leCount && (lidData <= eidData || justPrintFlag == 1))
        % -- print the landmark edge
        fprintf(outFile, 'EDGE_SE2_XY %d %d', landmarkEdgeData(1,lEIdx), landmarkEdgeData(2,lEIdx));
        for i = 3:7
            fprintf(outFile, ' %f', landmarkEdgeData(i, lEIdx));
        end
        fprintf(outFile,'\n');
        % -- increment the landmark edge index to read the next one
        lEIdx = lEIdx + 1;
    elseif ( eIdx <= eCount && (lidData >= eidData || justPrintFlag == 1))
        % -- print the edge
        fprintf(outFile, 'EDGE_SE2 %d %d', edgeData(1,eIdx), edgeData(2, eIdx));
        for i = 3:11
            fprintf(outFile,' %f', edgeData(i,eIdx));
        end
        fprintf(outFile,'\n');
        % -- increment the edge index to read the next one
        eIdx = eIdx + 1;
    end
    if (lEIdx > leCount && eIdx > eCount)
        break;
    end
end
fclose(outFile);
        
end