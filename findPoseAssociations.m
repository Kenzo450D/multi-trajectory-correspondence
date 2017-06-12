function [poseAssc, gtDenseAssc] = findPoseAssociations(fileName1, fileName2)
%FINDLANDMARKASSOCIATIAONS Find the landmark associations between two .mat
%files.
%  It loads the landmark data of both the files. Compares which ones are
%  close together, and maps them in a matrix.
% -------------------------------------------------------------------------
%  Input:
%    fileName1: filename including path for the mat file
%    fileName2: filename including path for the mat file (containing second
%    part of information)
% -------------------------------------------------------------------------
%  Output:
%    landmarkAssc: Associations between landmarks. It is a vector of the
%    size of that of the first landmark count. It is initialized to zero,
%    and only the values which correspond to the landmark of the other file
%    is declared as non-zero value (equals to the landmark number in the
%    other file)
% -------------------------------------------------------------------------

%% load the files
% -- load file 1
load(fileName1);
vertices1 = vertices;
landmarks1 = landmarks;
lCount1 = lCount;
vCount1 = vCount;
% -- load file 2
load(fileName2);
vertices2 = vertices;
landmarks2 = landmarks;
lCount2 = lCount;
vCount2 = vCount;

% total vertex and index adjustment
totalVertexCount1 = lCount1 + vCount1;

%% loop through the first landmarks
poseAssc = zeros(1,vCount1);
gtDenseAssc = zeros(1,totalVertexCount1);
for i = 1:vCount1
    vx = vertices1(i).x;
    vy = vertices1(i).y;
    for j = 1:vCount2
        vx2 = vertices2(j).x;
        vy2 = vertices2(j).y;
        diffx = vx - vx2;
        diffy = vy - vy2;
        if ((diffx ==0)  && (diffy == 0))
            poseAssc(i) = j;
            gtDenseAssc(vertices1(i).id) = vertices2(j).id;
        end
    end
end

end