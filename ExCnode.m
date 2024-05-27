function [x,y,z] = ExCnode(index,xyz)
% input element nodes connet info
% output the coordinate of nodes
len = length(index);
x = zeros(len,1);
y = zeros(len,1);
z = zeros(len,1);

for counti = 1:length(index)
    if index(counti)~=xyz(index(counti),1)
        fprintf(2,"Node Number may NOT start from 1!\n");
    end
    x(counti) = xyz(index(counti),2);
    y(counti) = xyz(index(counti),3);
    if length(xyz(1,:))>3
        z(counti) = xyz(index(counti),4);
    end
end
end