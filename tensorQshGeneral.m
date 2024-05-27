function [Qsh,Qsh6]=tensorQshGeneral(Jr,Js,Jt,Lx,Ly,Lz)

format long;
% Xr = Jr(1)*Lx(1)+Jr(2)*Lx(2)+Jr(3)*Lx(3);
% Yr = Jr(1)*Ly(1)+Jr(2)*Ly(2)+Jr(3)*Ly(3);
% Zr = Jr(1)*Lz(1)+Jr(2)*Lz(2)+Jr(3)*Lz(3);
% Xs = Js(1)*Lx(1)+Js(2)*Lx(2)+Js(3)*Lx(3);
% Ys = Js(1)*Ly(1)+Js(2)*Ly(2)+Js(3)*Ly(3);
% Zs = Js(1)*Lz(1)+Js(2)*Lz(2)+Js(3)*Lz(3);
% Xt = Jt(1)*Lx(1)+Jt(2)*Lx(2)+Jt(3)*Lx(3);
% Yt = Jt(1)*Ly(1)+Jt(2)*Ly(2)+Jt(3)*Ly(3);
% Zt = Jt(1)*Lz(1)+Jt(2)*Lz(2)+Jt(3)*Lz(3);

Xr = dot(Jr,Lx);
Yr = dot(Jr,Ly);
Zr = dot(Jr,Lz);
Xs = dot(Js,Lx);
Ys = dot(Js,Ly);
Zs = dot(Js,Lz);
Xt = dot(Jt,Lx);
Yt = dot(Jt,Ly);
Zt = dot(Jt,Lz);

% E_{ij}*e^{i}*e^{j} = vE_{kl}*j^{k}*j^{l};
% E_{ij} = Qsh_{ikjl}*vE_{kl};

Qsh = [Xr*Xr,Xs*Xs,Xt*Xt,Xs*Xt,Xt*Xs,Xr*Xt,Xt*Xr,Xr*Xs,Xs*Xr;
       Yr*Yr,Ys*Ys,Yt*Yt,Ys*Yt,Yt*Ys,Yr*Yt,Yt*Yr,Yr*Ys,Ys*Yr;
       Zr*Zr,Zs*Zs,Zt*Zt,Zs*Zt,Zt*Zs,Zr*Zt,Zt*Zr,Zr*Zs,Zs*Zr;
       Yr*Zr,Ys*Zs,Yt*Zt,Ys*Zt,Yt*Zs,Yr*Zt,Yt*Zr,Yr*Zs,Ys*Zr;
       Zr*Yr,Zs*Ys,Zt*Yt,Zs*Yt,Zt*Ys,Zr*Yt,Zt*Yr,Zr*Ys,Zs*Yr;
       Xr*Zr,Xs*Zs,Xt*Zt,Xs*Zt,Xt*Zs,Xr*Zt,Xt*Zr,Xr*Zs,Xs*Zr;
       Zr*Xr,Zs*Xs,Zt*Xt,Zs*Xt,Zt*Xs,Zr*Xt,Zt*Xr,Zr*Xs,Zs*Xr;
       Xr*Yr,Xs*Ys,Xt*Yt,Xs*Yt,Xt*Ys,Xr*Yt,Xt*Yr,Xr*Ys,Xs*Yr;
       Yr*Xr,Ys*Xs,Yt*Xt,Ys*Xt,Yt*Xs,Yr*Xt,Yt*Xr,Yr*Xs,Ys*Xr;];

% 转换成6X6的形式
proj = Qsh(1:3,:);
proj(4,:) = Qsh(4,:)+ Qsh(5,:);
proj(5,:) = Qsh(6,:)+ Qsh(7,:);
proj(6,:) = Qsh(8,:)+ Qsh(9,:);
Qsh6 = proj(:,1:3);
Qsh6(:,4) = proj(:,4)+proj(:,5);
Qsh6(:,5) = proj(:,6)+proj(:,7);
Qsh6(:,6) = proj(:,8)+proj(:,9);
Qsh6(:,4:6)=Qsh6(:,4:6)/2.0;

end