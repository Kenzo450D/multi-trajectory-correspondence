function [ heatDistance ] = distanceEmbedding( eigenVectors,eigenValues,nEigenValues, t_scale)

%% Heat Embedding
% figure();
% plotData = exp(-(t_scale/2));
% plot(plotData,eigenValues);

Y=eigenVectors(:,1:nEigenValues) * diag(exp(-(t_scale/2)*eigenValues(1:nEigenValues)));
heatDistance = Y*Y';

end