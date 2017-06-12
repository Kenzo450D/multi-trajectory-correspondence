function [vd] = amplifyHeatDescriptor(poseNeighbourhoodMap, vCount, vd)
%AMPLIFYHEATDESCRIPTOR Amplifies a vertex descriptor based on the
%poseNeighbourhood Map.

%% parameter initialization
increaseParam1 = 10; % multiplier for the tScale chosen
increaseParam2 = 3; % multiplier for the neighbours (1 hop) of the tScale chosen
increaseParam3 = 2;
baseParam = 0.9; % multiplier for the rest of the tScales
ntScales = size(vd,1);

%% get the neighbourhood type
[poseNeighbourhoodType] = declareNeighbourhoodSparseDense(poseNeighbourhoodMap, vCount);

%% use the neighbourhood type to amplify the pose descriptors

amplifier = ones(ntScales,1)*baseParam;
amplifier(ntScales) = increaseParam1;
amplifier(ntScales-1) = increaseParam2;
amplifier(ntScales-2) = increaseParam3;
matAmp = repmat(amplifier',size(vd,2),1)';
for i = 1:vCount
    if(poseNeighbourhoodType(i) == 1)
        % -- sparse locality... needs to be amplified
        vd(:,:,i) = vd(:,:,i).*matAmp;
    end
end


end