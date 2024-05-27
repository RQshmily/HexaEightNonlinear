function []=mplot(Value,gcoords,coords,nodes)
map=[0,0,255;
    0,93,255;
    0,185,255;

    0,255,232;
    0,255,139;

    0,255,46;
    46,255,0;

    139,255,0;
    232,255,0;

    255,185,0;
    255,93,0;
    255,0,0];

nes = length(nodes(:,1));
figure(66);
clf;
hold on;
[faces,parr]=getfaces();


for count=1:nes
    cnode = nodes(count,:);


    [x0,y0,z0] = ExCnode(cnode,gcoords);

    xp0=x0(parr);
    yp0=y0(parr);
    zp0=z0(parr);

    plot3(xp0,yp0,zp0,'--k');
    hold on;


    [xcoord,ycoord,zcoord] = ExCnode(cnode,coords);

    Pd=Value(cnode);

    for c=1:length(faces)
        subf=faces{c};

        subx=xcoord(subf);
        suby=ycoord(subf);
        subz=zcoord(subf);
        subv=Pd(subf);

        patch(subx,suby,subz,subv);
    end
end



axis tight;
colormap(map/255);
axis equal
set(gca,'fontsize',12);
% colorbar;
minv=min(Value);
maxv=max(Value);
if minv==maxv
    minv=minv-0.5;
    maxv=maxv+0.5;
end

caxis([minv,maxv]);% 设置颜色数据为最大值和最小值
t1=caxis;
t1=linspace(t1(1),t1(2),12);
colorbar('ytick',t1);
view(3);

end

function [faces,parr]=getfaces()

subf1=[1,2,3,4];
subf2=[5,6,7,8];
subf3=[1,2,6,5];
subf4=[2,3,7,6];
subf5=[3,4,8,7];
subf6=[1,4,8,5];
parr = [1,2,3,4,1,5,6,2,6,7,3,7,8,4,8,5];
faces={subf1;subf2;subf3;subf4;subf5;subf6;};

end