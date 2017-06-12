function [eucDiff] = getPoseDiffCombinations(vData)
%GETPOSEDIFFCOMBINATIONS Get euclidean difference of all combinations.

vCount = size(vData,2);
a = [1:vCount];
xData = zeros(vCount);
yData = zeros(vCount);
for i = 1:vCount
    for j = i+1:vCount
        xData(i,j) = vData(2,i) - vData(2,j);
        yData(i,j) = vData(3,i) - vData(3,j);        
    end
end
xData = xData + xData';
yData = yData + yData';
xData = xData.^2;
yData = yData.^2;
eucDiff = xData + yData;
eucDiff = eucDiff.^0.5;
end