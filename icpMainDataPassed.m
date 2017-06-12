function icpMainDataPassed(vertices1, vertices2, gtPoseAssc, atePrecisionFile, outMatICPFileName)
%ICPMAINDATAPASSED It is a data passed version of ICP Only code. This
%accepts the vertex and edge information, and runs unsupervised ICP on it. 

%% call the ICP function

[Ricp, Ticp, ER, t, newMatchesMade] = runICPOnly(vertices1, vertices2);

% -- convert to matrix
vData1 = getVertexMatrixForm(vertices1);
vData2 = getVertexMatrixForm(vertices2);

%% ate and precision

% -- convert gtPoseAssc to gtMatchesMade
poseIdx = gtPoseAssc>0;
posesToCompare = gtPoseAssc(poseIdx);
gtMatchesMade = zeros(length(posesToCompare),2);
for i = 1:length(posesToCompare)
    if (poseIdx(i) ~= 0)
        gtMatchesMade(i,1) = i;
        gtMatchesMade(i,2) = gtPoseAssc(i);
    end
end

[ate] = calculateATE(vData1, vData2, newMatchesMade);
[p, relaxPrecision] = calculatePrecision(gtMatchesMade, newMatchesMade, (size(vData1,2)*0.01));

outFile = atePrecisionFile;
fid = fopen(outFile,'a')';

fprintf(1,'ICP: %f %f %f\n',ate,p,relaxPrecision);
fprintf(fid,'%f %f %f\n',ate,p,relaxPrecision);
fclose(fid);

%% save the matches made in the file
icpMatchesMade = newMatchesMade;
save(outMatICPFileName,'icpMatchesMade');

end