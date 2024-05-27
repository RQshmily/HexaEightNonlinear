%==========================================================================%
%                  HEAD                                                    %
%               Implicit_Nonlinear                                         %
%          RANQI 2023-10-27                                                %
%              HUNAN UNIVERSITY.                                           %
%                  ranqi@hnu.edu.cn                                        %
%==========================================================================%
%                                                                          %
%-------------------------- Geometric nonlinear ---------------------------%


clear;
close all;
clc;
format long;

load("HS8_Solid.mat");

%==========================================================================

% PAPERING WORK
nstep = 25;
CurCoords = gcoords;
tol = 1.0e-4;% 两个迭代步之间的位移变化量
Set_ZERO = zeros(3*length(gcoords(:,1)),1);
totaldisp = Set_ZERO;
dispinr = Set_ZERO;

fext = Set_ZERO;

% inputf size:n X 2
% inputf(:,1) load applied dofs
% inputf(:,2) corresponding load values
fext(inputf(:,1))=inputf(:,2);

nes = length(nodes(:,1));
nns = length(gcoords(:,1));             % number of nodes of system.


% stress--updated pk2
for co = 1:nes
    for ip = 1:8
        info.stress{co,ip} = zeros(6,1);
    end
end

Cm = fematiso(4,Mat{2,2},Mat{3,2});
nds = 3 * nns;                       % number of DOFs of system.

totSteps = 0;
disp = Set_ZERO;

for istep = 1:nstep
    factor=istep/nstep;

    fprintf("======= %3d =======\n",factor);


    dispinr = Set_ZERO;


    conv = 0;iter = 0;
    kk = zeros(nds,nds);
    while ~conv

        [ kk,fint,info ] = getKF_C3D8_Bathe(gcoords,nodes,Cm,dispinr,disp,info);



        kk(BC(:,1),:)=0;
        kk(:,BC(:,1))=0;
        for ii=1:length(BC(:,1))
            kk(BC(ii,1),BC(ii,1))=1.0;
        end

        rf = factor*fext-fint;

        rf(BC(:,1)) = factor*BC(:,2)-disp(BC(:,1));

        dispinr = sparse(kk)\rf;
        disp = disp+dispinr;

        Inr=max(abs(dispinr));
        mrf=max(abs(rf));
        if iter==0
            tolf = mrf/1000;
            Tar = Inr/10000;
        end

        fprintf("%3d  Rf: %2.3e  Disp:%2.3e \n",iter,mrf,Inr);
        iter=iter+1;
        if (mrf<tolf)&&(Inr<Tar)
            conv = 1;
        end
    end
    totSteps = totSteps+iter;

    % 更新当前坐标
    CurCoords(:,2) = gcoords(:,2)+disp(1:3:end);
    CurCoords(:,3) = gcoords(:,3)+disp(2:3:end);
    CurCoords(:,4) = gcoords(:,4)+disp(3:3:end);

    x = disp(1:ndpn:end);
    y = disp(2:ndpn:end);
    z = disp(3:ndpn:end);
    M = sqrt(x.*x + y.*y + z.*z);
    mplot(M,gcoords,CurCoords,nodes);
    pause(0.001);

end

% End
fprintf(2,"===================== Analysis finished =====================\n");
fprintf(2,"===================== Total Steps: %d =====================\n",totSteps);

