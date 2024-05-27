function [ stiff,fint,info ] = getKF_C3D8(gcoords,nodes,D,dispinr,totdisp,info)
format long;
nnpe = 8;
ndpn = 3;
nes = length(nodes(:,1));
nns = length(gcoords(:,1));

stiff = zeros(3*nns,3*nns);
fint = zeros(3*nns,1);

point(1) = -sqrt(3)/3.0;
point(2) = -point(1);
dnds = zeros(24,8);
lip = 0;
for countx=1:2
    ksi=point(countx);
    for county=1:2
        eta=point(county);
        for countz=1:2
            zeta=point(countz);
            lip = lip+1;
            [dnds_ip,~]=getdNds(ksi,eta,zeta,'Hexa1');
            dnds(3*lip-2:3*lip,:) = dnds_ip;
        end
    end
end
I3 = eye(3);



for iel=1:nes
    ke = zeros(24,24);
    f = zeros(24,1);

    cnode = nodes(iel,:);

    index = feeldof(cnode,nnpe,ndpn);
    edispt = [totdisp(3*cnode-2),totdisp(3*cnode-1),totdisp(3*cnode-0)];
    edispt([2,3,6,7],1)=1.0;

    cord = gcoords(cnode,2:4);
    for ip = 1:8
        Jac = dnds(3*ip-2:3*ip,:)*cord;
        dN_dX = Jac\dnds(3*ip-2:3*ip,:);
        F = dN_dX*edispt + I3;

        for inode = 1:8
            B_L(:,(inode*3-2:inode*3)) = [F(1,1)*dN_dX(1,inode),F(2,1)*dN_dX(1,inode),F(3,1)*dN_dX(1,inode);
                F(1,2)*dN_dX(2,inode),F(2,2)*dN_dX(2,inode),F(3,2)*dN_dX(2,inode);
                F(1,3)*dN_dX(3,inode),F(2,3)*dN_dX(3,inode),F(3,3)*dN_dX(3,inode);
                F(1,1)*dN_dX(2,inode)+F(1,2)*dN_dX(1,inode),F(2,1)*dN_dX(2,inode)+F(2,2)*dN_dX(1,inode),F(3,1)*dN_dX(2,inode)+F(3,2)*dN_dX(1,inode);
                F(1,2)*dN_dX(3,inode)+F(1,3)*dN_dX(2,inode),F(2,2)*dN_dX(3,inode)+F(2,3)*dN_dX(2,inode),F(3,2)*dN_dX(3,inode)+F(3,3)*dN_dX(2,inode);
                F(1,3)*dN_dX(1,inode)+F(1,1)*dN_dX(3,inode),F(2,3)*dN_dX(1,inode)+F(2,1)*dN_dX(3,inode),F(3,3)*dN_dX(1,inode)+F(3,1)*dN_dX(3,inode)];

            B_NL(:,(inode*3-2:inode*3))= [dN_dX(1,inode),0,0;
                dN_dX(2,inode),0,0;
                dN_dX(3,inode),0,0;
                0,dN_dX(1,inode),0;
                0,dN_dX(2,inode),0;
                0,dN_dX(3,inode),0;
                0,0,dN_dX(1,inode);
                0,0,dN_dX(2,inode);
                0,0,dN_dX(3,inode);];
        end

%         R = 0.5*(F'*F-I3);
%         S6 = D*m_2_v6(R);
        L = 0.5*(F*F'-I3);
        S6 = D*m_2_v6(L);
        S3 = v6_2_m(S6);

        I = blkdiag(S3,S3,S3);

        ke  = ke + det(Jac)*(B_L'*D*B_L+B_NL'*I*B_NL);
        f = f + det(Jac)*(B_L'*[S3(1,1);S3(2,2);S3(3,3);S3(1,2);S3(2,3);S3(1,3)]);
    end


    index = feeldof(cnode,nnpe,3);
    stiff = assemblekk(stiff,ke,index,index);
    fint = assemblekk(fint,f,index,1);

end

% for idsi = 1:nnpe
% 
%     for idsj = 1:nnpe
%         stiff(nodes(:,idsi)*3-2,nodes(:,idsj)*3-2) = stiff(nodes(:,idsi)*3-2,nodes(:,idsj)*3-2) + ke(:,idsi*3-2,idsj*3-2);
%         stiff(nodes(:,idsi)*3-2,nodes(:,idsj)*3-1) = stiff(nodes(:,idsi)*3-2,nodes(:,idsj)*3-1) + ke(:,idsi*3-2,idsj*3-1);
%         stiff(nodes(:,idsi)*3-2,nodes(:,idsj)*3-0) = stiff(nodes(:,idsi)*3-2,nodes(:,idsj)*3-0) + ke(:,idsi*3-2,idsj*3-0);
%         stiff(nodes(:,idsi)*3-1,nodes(:,idsj)*3-2) = stiff(nodes(:,idsi)*3-1,nodes(:,idsj)*3-2) + ke(:,idsi*3-1,idsj*3-2);
%         stiff(nodes(:,idsi)*3-1,nodes(:,idsj)*3-1) = stiff(nodes(:,idsi)*3-1,nodes(:,idsj)*3-1) + ke(:,idsi*3-1,idsj*3-1);
%         stiff(nodes(:,idsi)*3-1,nodes(:,idsj)*3-0) = stiff(nodes(:,idsi)*3-1,nodes(:,idsj)*3-0) + ke(:,idsi*3-1,idsj*3-0);
%         stiff(nodes(:,idsi)*3-0,nodes(:,idsj)*3-2) = stiff(nodes(:,idsi)*3-0,nodes(:,idsj)*3-2) + ke(:,idsi*3-0,idsj*3-2);
%         stiff(nodes(:,idsi)*3-0,nodes(:,idsj)*3-1) = stiff(nodes(:,idsi)*3-0,nodes(:,idsj)*3-1) + ke(:,idsi*3-0,idsj*3-1);
%         stiff(nodes(:,idsi)*3-0,nodes(:,idsj)*3-0) = stiff(nodes(:,idsi)*3-0,nodes(:,idsj)*3-0) + ke(:,idsi*3-0,idsj*3-0);
% 
%     end
%     fint(nodes(:,idsi)*3-2) = fint(nodes(:,idsi)*3-2) + f(:,idsi*3-2);
%     fint(nodes(:,idsi)*3-1) = fint(nodes(:,idsi)*3-1) + f(:,idsi*3-1);
%     fint(nodes(:,idsi)*3-0) = fint(nodes(:,idsi)*3-0) + f(:,idsi*3-0);
% end



end

function [v] = m_2_v9(m)

v=zeros(9,1);
v(1)=m(1,1); v(4)=m(1,2); v(6)=m(1,3);
v(7)=m(2,1); v(2)=m(2,2); v(5)=m(2,3);
v(9)=m(3,1); v(8)=m(3,2); v(3)=m(3,3);

end

function [m] = v9_2_m(v)

m = [v(1),v(4),v(6);
     v(7),v(2),v(5);
     v(9),v(8),v(3);];

end

function [v] = m_2_v6(m)

v = zeros(6,1);
v(1) = m(1,1);
v(2) = m(2,2);
v(3) = m(3,3);
v(4) = (m(1,2) + m(2,1));
v(5) = (m(2,3) + m(3,2));
v(6) = (m(1,3) + m(3,1));


% v(4) = (m(1,2) + m(2,1))/2;
% v(5) = (m(2,3) + m(3,2))/2;
% v(6) = (m(1,3) + m(3,1))/2;

end

function [m] = v6_2_m(v)

m = [v(1),v(4)/2,v(6)/2;
     v(4)/2,v(2),v(5)/2;
     v(6)/2,v(5)/2,v(3);];

% m = [v(1),v(4),v(6);
%      v(4),v(2),v(5);
%      v(6),v(5),v(3);];

end
