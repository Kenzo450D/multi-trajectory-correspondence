function [stepDiff] = getPoseStepDiffCombinations(vData)
%GETPOSEDIFFCOMBINATIONS Get step difference of all combinations.

vCount = size(vData,2);
stepDiff = zeros(vCount, vCount);
for i = 1:vCount
    for j = i+1:vCount
        stepDiff(i,j) = abs(j - i);
    end
end
stepDiff = stepDiff + stepDiff';
end