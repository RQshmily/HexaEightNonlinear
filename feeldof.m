function [index] = feeldof(node,nnpe,ndpn)
% �ҵ���Ԫ�ڵ�node����Ӧ�����ɶȵ�λ��
% node ��Ԫ�ڵ���Ϣ
% nnpe ��Ԫ�ڵ���
% ndpn �ڵ����ɶ���
index =zeros(nnpe*ndpn,1);
for count = 1:nnpe
    for countin = 1:ndpn
        index(count*ndpn-ndpn+countin) = ndpn*node(count)-(ndpn-countin);
    end
end
end