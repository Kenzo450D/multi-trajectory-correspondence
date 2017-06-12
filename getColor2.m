function [color2] = getColor2(color1, vCount2, denseAssc)

color2 = zeros(vCount2, 3);
color2(:,2) = 1;
for i = 1:length(denseAssc)
    if denseAssc(i) ~= 0
        color2(denseAssc(i),:) = color1(i,:);
    end
end
end