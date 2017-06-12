function [R, T] = icpForMatchesMade(matchesMade, vData1, vData2)



%% convert matchesMade to the point cloud as is required by the icp code
matchCount = size(matchesMade,1);
% ---- convert vData to just the required two information
poseData1 = vData1(2:3,:);
poseData2 = vData2(2:3,:);
icpInputData1 = zeros(3,matchCount);
icpInputData2 = zeros(3,matchCount);
% matchesMade is a 2column matrix.

for i = 1:matchCount
    % -- read each match at a time, and store the information
    icpInputData1(1:2,i) = poseData1(:,matchesMade(i,1));
    icpInputData2(1:2,i) = poseData2(:,matchesMade(i,2));
end
q = icpInputData1;
p = icpInputData2;

%% get the information or RT for the transformation
p_idx = true(1,length(p));
weights = ones(1,size(matchesMade, 1));
[R,T] = eq_point(q,p, weights(p_idx));

%% Apply last transformation
%pt = R * p + repmat(TT(:,:,k+1), 1, Np);

end

function [R,T] = eq_point(q,p,weights)

m = size(p,2);
n = size(q,2);

% normalize weights
weights = weights ./ sum(weights);

% find data centroid and deviations from centroid
q_bar = q * transpose(weights);
q_mark = q - repmat(q_bar, 1, n);
% Apply weights
q_mark = q_mark .* repmat(weights, 3, 1);

% find data centroid and deviations from centroid
p_bar = p * transpose(weights);
p_mark = p - repmat(p_bar, 1, m);
% Apply weights
%p_mark = p_mark .* repmat(weights, 3, 1);

N = p_mark*transpose(q_mark); % taking points of q in matched order

[U,~,V] = svd(N); % singular value decomposition

R = V*diag([1 1 det(U*V')])*transpose(U);

T = q_bar - R*p_bar;
end