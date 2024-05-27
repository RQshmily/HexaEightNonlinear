function [index] = feeldof(node,nnpe,ndpn)
% 找到单元节点node所对应的自由度的位置
% node 单元节点信息
% nnpe 单元节点数
% ndpn 节点自由度数
index =zeros(nnpe*ndpn,1);
for count = 1:nnpe
    for countin = 1:ndpn
        index(count*ndpn-ndpn+countin) = ndpn*node(count)-(ndpn-countin);
    end
end
end