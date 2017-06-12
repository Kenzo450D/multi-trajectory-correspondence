function [agVertexList] = getAssociationsVertexSubset(vd1, vd2, v1start, v1end, v2start, v2end, neighbourLimit, totalVertexCount1, tScale1, tScale2, poseNeighbourhoodMap1, poseNeighbourhoodMap2)
% -- Input
%   vd1 : vertex descriptor set 1
%   vd2 : vertex descriptor set 2
%   v1start : start boundary of vertex set 1
%   v1end : end boundary of vertex set 1
%   v2start : start boundary of vertex set 2
%   v2end : end boundary of vertex set 2
%   neighbourLimit: number of neighbours to be considered in Association Graph
%   totalVertexCount1: total vertex count for the first graph, needed to arrange data
%   tScale1 : tScale parameters for graph 1
%   tScale2 : tScale parameters for graph 2
%   *POSENEIGHBOURHOODMAP* : It is to denote how many neighbours does that
%   particular landmark territory have, it is so that a higher scale of
%   diffusion can be given higher importance in spare areas, and a lower
%   scale of diffusion can be given importance in dense areas.
%   poseNeighbourhoodMap1 : pose Neighbourhood information for graph 1
%   poseNeighbourhoodMap2 : pose Neighbourhood information for graph 2)
%   poseNeighbourhoodMap1


% ---- make sparse matrix to make their adjacency graph

agVertexList = [];

% -- form the vertices
%tic;
descr_threshold = 0.05;
for i = v1start:v1end
    index1 = i;
    hdv1   = vd1(:,:,index1);
    agVertexTmpList = [];
    agIndexTmpList = [];
    cmpHds_arr = [];
    for j = v2start:v2end
        index2 = j;
        hdv2   = vd2(:,:,index2);
        % -- compare head descriptors if only they hold significant
        % information. And most of their values are not zeros.
        if (norm(hdv1) < 1e-6 || norm(hdv2) < 1e-6)
            continue;
        end
        % -- compare heat descriptors based on the difference of the
        % sections
        
        %cmpHd = getComparedHeatDescriptor(hdv1, hdv2, tScale1, tScale2, poseNeighbourhoodMap1(i), poseNeighbourhoodMap2(j));
        cmpHd = norm(hdv1 - hdv2); % compare between two heat descriptor for vertices
        if (cmpHd < descr_threshold)
            cmpHds_arr = [cmpHds_arr; cmpHd];
            tmpVal = [index1, index2];
            agVertexTmpList = [agVertexTmpList;tmpVal];
        end
    end
    [~, tindx1]=sort(cmpHds_arr);
     if length(cmpHds_arr)>neighbourLimit
        for j = 1:neighbourLimit
            agVertexList = [agVertexList; agVertexTmpList(tindx1(j),:)];
        end
     elseif length(cmpHds_arr)>0
            agVertexList = [agVertexList; agVertexTmpList];
     %else
     %       disp(sprintf('No Match Found for index: %d in left trajectory',index1));
     end
    %fprintf(1,'length of agVertexList: %d\n',length(agVertexList));
end

agVertexCount = length(agVertexList);
%fprintf(1,'Time taken to create vertices: %d\n',toc);

end