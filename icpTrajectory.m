function [ match, mindist ] = icpTrajectory( q, p, initMatch )
%ICPTRAJECTORY Trajectory 2 is matched with Trajectory 1
%   Brute Force matching between the two trajectories. It is considered
%   that a few of the 
    m = size(p,2);
    n = size(q,2);  
    if (exits('initMatch'))
        match = initMatch;
    else
        match = zeros(1,m);
    end
    mindist = zeros(1,m);
    for ki=1:m
        d=zeros(1,n);
        for ti=1:2
            d=d+(q(ti,:)-p(ti,ki)).^2;
        end
        [mindist(ki),match(ki)]=min(d);
    end
    mindist = sqrt(mindist);
end

% Code copied from
% function [match mindist] = match_bruteForce(q, p)
%     m = size(p,2);
%     n = size(q,2);    
%     match = zeros(1,m);
%     mindist = zeros(1,m);
%     for ki=1:m
%         d=zeros(1,n);
%         for ti=1:3
%             d=d+(q(ti,:)-p(ti,ki)).^2;
%         end
%         [mindist(ki),match(ki)]=min(d);
%     end
%     
%     mindist = sqrt(mindist);

