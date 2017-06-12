function [result] = quantitativeTest(gtAssc, poseAssc)
%QUANTITATIVETEST A Quantitative analysis of quality of results

result = zeros(size(gtAssc));
result = result + 2;
numberAssc = nnz(gtAssc); % number of total associations
for i = 1:length(gtAssc)
    if (gtAssc(i) ~= 0)
        result(i) = (abs(poseAssc(i) - gtAssc(i)))/numberAssc;
    elseif (poseAssc(i) == 0)
        result(i) = 3;
    end
end
end