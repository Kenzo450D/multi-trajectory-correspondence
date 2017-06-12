function convertToMatFile(vData, lData, eData, leData, outFile, outG2oFileName)

%% convert the data to structure

% -- get the count
vCount = size(vData,2);
eCount = size(eData, 2);
lCount = size(lData,2);
leCount = size(leData,2);

for i = 1:vCount
    vertices(i).id=vData(1,i);
    vertices(i).x=vData(2,i);
    vertices(i).y=vData(3,i);
    vertices(i).o=vData(4,i);
end
for i = 1:lCount
    landmarks(i).id=lData(1,i);
    landmarks(i).x=lData(2,i);
    landmarks(i).y=lData(3,i);
end

for i = 1:eCount
    edges(i).v1   = eData(1,i);
    edges(i).v2   = eData(2,i);
    edges(i).dx   = eData(3,i);
    edges(i).dy   = eData(4,i);
    edges(i).dth  = eData(5,i);
    %fill Covariance Matrix
    covMatrix          = zeros(3);
    covMatrix(1,1)     = eData(6,i);
    covMatrix(1,2)     = eData(7,i);
    covMatrix(2,1)     = eData(7,i);
    covMatrix(1,3)     = eData(8,i);
    covMatrix(3,1)     = eData(8,i);
    covMatrix(2,2)     = eData(9,i);
    covMatrix(2,3)     = eData(10,i);
    covMatrix(3,2)     = eData(10,i);
    covMatrix(3,3)     = eData(11,i);
    edges(i).covMatrix = covMatrix;
end
for i = 1:leCount
    ledges(i).v1 = leData(1,i);
    ledges(i).v2 = leData(2,i);
    ledges(i).dx = leData(3,i);
    ledges(i).dy = leData(4,i);
    % fill covariance matrix
    covMatrix  = zeros(2);
    covMatrix(1,1) = leData(5,i);
    covMatrix(1,2) = leData(6,i);
    covMatrix(2,1) = leData(6,i);
    covMatrix(2,2) = leData(7,i);
    ledges(i).covMatrix = covMatrix;
end
if (~isempty(outG2oFileName))
    exportLandmarkG2oFile(outG2oFileName, vData, lData, eData, leData);
    %exportLandmarkG2oGraph(vertices,landmarks, edges, ledges, outG2oFileName);
end
dim = 2;
save(outFile,'vertices','vCount','landmarks', 'lCount','edges','eCount','ledges','leCount', 'dim');

end