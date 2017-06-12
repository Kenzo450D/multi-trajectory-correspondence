function [p, relaxPrecision] = calculatePrecision(gtMatchesMade, autMatchesMade, thresholdRelax)
%CALCULATEPRECISION Calculates precision in matching
% Input: gtMatchesMade and autMatchesMade (2 column array)
% Output:
%    p: precision
%    relaxPrecision: relaxed precision in which even if the index is in a
%    neighbourhood of the original pose, it would still work.

%precision = truePos / (truePos + falsePos)

% -- initialization
%thresholdRelax = 5;


% -- arrange autMatches sorted according to the first column
[~,sortIdx] = sort(autMatchesMade(:,1));
autMatchesMade(:,1) = autMatchesMade(sortIdx,1);
autMatchesMade(:,2) = autMatchesMade(sortIdx,2);

% -- make the true positive
% truePos are the ones which would have the same index of vData2, as
% gtPoseAssc

tpCount = 0;
fpCount = 0;
rlxTpCount = 0;
rlxFpCount = 0;

for i = 1:length(gtMatchesMade)
    gtv1Val = gtMatchesMade(i,1);
    gtv2Val = gtMatchesMade(i,2);
    % -- get vertex for vData2
    autv1Idx = find(autMatchesMade(:,1) == gtv1Val);
    if (~isempty(autv1Idx))
        l = length(autv1Idx);
        for j = 1:l
            autv2Val = autMatchesMade(autv1Idx(j),2);
            if (autv2Val == gtv2Val)
                tpCount = tpCount + 1;
            else
                fpCount = fpCount + 1;
            end
            diffn = abs(gtv2Val - autv2Val);
            if (diffn <= thresholdRelax)
                rlxTpCount = rlxTpCount + 1;
            else
                rlxFpCount = rlxFpCount + 1;
            end
        end
    end
end

p = tpCount / (tpCount + fpCount);
relaxPrecision = rlxTpCount / (rlxTpCount + rlxFpCount);
end


% % -- make true positive with relaxation of threshold steps
% 
% rlxTpCount = 0;
% rlxFpCount = 0;
% 
% for i = 1:length(gtMatchesMade)
%     gtv1Val = gtMatchesMade(i,1);
%     gtv2Val = gtMatchesMade(i,2);
%     % -- get vertex for vData2
%     autv1Idx = find(autMatchesMade(:,1) == gtv1Val);
%     if (~isempty(autv1Idx))
%         l = length(autv1Idx);
%         for j = 1:l
%             autv2Val = autMatchesMade(autv1Idx(j),2);
%             diffn = abs(gtv2Val - autv2Val);
%             fprintf(1,'diffn: %d\n',diffn);
%             if (diffn <= thresholdRelax)
%                 rlxTpCount = rlxTpCount + 1;
%             else
%                 rlxFpCount = rlxFpCount + 1;
%             end
%         end
%     end
% end
% 
% relaxPrecision = rlxTpCount / (rlxTpCount + rlxFpCount);
% end