function [kk]= assemblekk(kk,k,indexi,indexj)
% 该函数进行单元刚度矩阵的组装
% 输入单元刚度矩阵 k
% 输入单元的自由度index
% 输出整体刚度矩阵 KK
%---------------------------------------------------------------
Leni = length(indexi);
Lenj = length(indexj);

for n1=1:Leni
    for n2=1:Lenj
        kk(indexi(n1),indexj(n2))= kk(indexi(n1),indexj(n2))+k(n1,n2);
    end
end