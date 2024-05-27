function [ stif,fint,info ] = getKF_C3D8_Bathe(gcoords,nodes,D,dispinr,totdisp,info)
format long;
% Calculate strain at covariant basis
% Incremental TL
% ranqi@hnu.edu.cn 
nnpe = 8;
ndpn = 3;
nes = length(nodes(:,1));
nns = length(gcoords(:,1));

stif = zeros(ndpn*nns,ndpn*nns);
fint = zeros(ndpn*nns,1);

point(1) = -sqrt(3)/3.0;
point(2) = -point(1);



for iel=1:nes
    ke = zeros(24,24);
    elf = zeros(24,1);

    cnode = nodes(iel,:);
    index = feeldof(cnode,nnpe,ndpn);
    edisptv = totdisp(index);
    edispiv = dispinr(index);

    edisptm = [edisptv(1:3:end),edisptv(2:3:end),edisptv(3:3:end)];
    edispim = [edispiv(1:3:end),edispiv(2:3:end),edispiv(3:3:end)];

    cord = gcoords(cnode,2:4);
    celec = cord + edisptm;

    lip = 0;
    for lpr=1:2
        r = point(lpr);
        for lps=1:2
            s = point(lps);
            for lpt=1:2
                t = point(lpt);

                lip = lip+1;
                [Blin,Bnon] = calBmat_C3D8_Bathe_UpdatedTL(r,s,t,celec,edisptm);

                dnds = getdNds(r,s,t);
                Jac = dnds*cord;
                ivJ = inv(Jac);
                iJ1 = ivJ(:,1);
                iJ2 = ivJ(:,2);
                iJ3 = ivJ(:,3);
                detJ0 = det(Jac);
                J1 = Jac(1,:);
                J2 = Jac(2,:);
                ez = cross(J1,J2);
                ex = cross(J2,ez);
                ey = cross(ez,ex);

                ex = ex/norm(ex);
                ey = ey/norm(ey);
                ez = ez/norm(ez);


                [~,Qsh]=tensorQshGeneral(iJ1,iJ2,iJ3,ex,ey,ez);


                strinr = Blin*edispiv;

                for ii=1:6
                    Bi = reshape(Bnon(ii,:,:),24,24);
                    strinr(ii) = strinr(ii) + 0.5*edispiv'*Bi*edispiv;
                end

                Cm = Qsh'*D*Qsh;

                steinr = Cm*strinr;

                ste = info.stress{iel,lip} + steinr;
                info.stress{iel,lip} = ste;

                ke = ke+ Blin' * Cm * Blin * detJ0;
                elf = elf + Blin' * ste * detJ0;

                for ii=1:6
                    Bi = reshape(Bnon(ii,:,:),24,24);
                    ke = ke + ste(ii).*Bi*detJ0;
                end

            end
        end
    end


    stif = assemblekk(stif,ke,index,index);
    fint = assemblekk(fint,elf,index,1);

end

end

