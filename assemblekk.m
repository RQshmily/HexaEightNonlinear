function [kk]= assemblekk(kk,k,indexi,indexj)
% �ú������е�Ԫ�նȾ������װ
% ���뵥Ԫ�նȾ��� k
% ���뵥Ԫ�����ɶ�index
% �������նȾ��� KK
%---------------------------------------------------------------
Leni = length(indexi);
Lenj = length(indexj);

for n1=1:Leni
    for n2=1:Lenj
        kk(indexi(n1),indexj(n2))= kk(indexi(n1),indexj(n2))+k(n1,n2);
    end
end