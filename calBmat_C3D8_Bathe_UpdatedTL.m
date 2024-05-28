function [Blin,Bnon] = calBmat_C3D8_Bathe_UpdatedTL(r,s,t,celec,U)
% 按照矢量完全形式，先求导后投影
format long;

[dndsAll]=getdNds(r,s,t);
dndr = dndsAll(1,:);
dnds = dndsAll(2,:);
dndt = dndsAll(3,:);

g1 = dndr*celec;
g2 = dnds*celec;
g3 = dndt*celec;

du1dr = zeros(3,24);
du1ds = zeros(3,24);
du1dt = zeros(3,24);

nnpe = 8;

for inode = 1:nnpe

    dnr = dndr(inode);
    dns = dnds(inode);
    dnt = dndt(inode);

    st = 3*(inode-1);
    du1dr(1,st + 1) =  dnr;
    du1dr(2,st + 2) =  dnr;
    du1dr(3,st + 3) =  dnr;

    du1ds(1,st + 1) =  dns;
    du1ds(2,st + 2) =  dns;
    du1ds(3,st + 3) =  dns;

    du1dt(1,st + 1) =  dnt;
    du1dt(2,st + 2) =  dnt;
    du1dt(3,st + 3) =  dnt;

end

e11 = g1*du1dr;
e22 = g2*du1ds;
e33 = g3*du1dt;
e12 = (g1*du1ds + g2*du1dr);
e13 = (g1*du1dt + g3*du1dr);
e23 = (g2*du1dt + g3*du1ds);
Blin = zeros(6,24);
Blin(1,:) = e11;
Blin(2,:) = e22;
Blin(3,:) = e33;
Blin(4,:) = e23;
Blin(5,:) = e13;
Blin(6,:) = e12;

Bnon = zeros(6,24,24);
Bnon(1,:,:) = du1dr'*du1dr;
Bnon(2,:,:) = du1ds'*du1ds;
Bnon(3,:,:) = du1dt'*du1dt;
Bnon(4,:,:) = du1ds'*du1dt + du1dt'*du1ds;
Bnon(5,:,:) = du1dr'*du1dt + du1dt'*du1dr;
Bnon(6,:,:) = du1dr'*du1ds + du1ds'*du1dr;

end

