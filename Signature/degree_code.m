function [code_set] = degree_code(com_fdeg,k1,c)

deg_code = cell2mat(com_fdeg(c,k1));
ck = size(k1,2);
code_book = ones(1,ck)';
for i=1:ck
    code_book(i) = 2^(i-1);
end
node_code = deg_code*code_book+1;
code_set = cell(2^ck,1);
for i=1:2^ck
    code_set{i} = find(node_code==i);
end
