function [edges] = removeLoopClosures(edges)
%REMOVELOOPCLOSURES Removes all edges in which the vertices are not
%successive. STRICT WARNING that is also removes landmark pose edges as
%well. It is advised that this function be ONLY used for loop closure
%graphs and NOT FOR LANDMARK POSE GRAPHS.

eCount = size(edges,2);
i = 1;
while(i <= eCount)
    v1 = edges(i).v1;
    v2 = edges(i).v2;
    diffn = abs(v2 - v1);
    if (diffn > 1)
        edges(i) = [];
        eCount = size(edges,2);
    else
        i = i + 1;
    end
end


end