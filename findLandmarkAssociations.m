function [landmarkAssc] = findLandmarkAssociations(fileName1, fileName2)
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
landmarks1 = landmarks;
lCount1 = lCount;
% -- load file 2
load(fileName2);
landmarks2 = landmarks;
lCount2 = lCount;

%% loop through the first landmarks
landmarkAssc = zeros(1,lCount1);
for i = 1:lCount1
    lx = landmarks1(i).x;
    ly = landmarks1(i).y;
    for j = 1:lCount2
        lx2 = landmarks2(j).x;
        ly2 = landmarks2(j).y;
        diffx = lx - lx2;
        diffy = ly - ly2;
        if ((diffx ==0)  && (diffy == 0))
            landmarkAssc(i) = j;
        end
    end
end

end