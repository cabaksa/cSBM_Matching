function [G,com_idx] = balanced_sbm(n,p,q,k)

% k = # of community
% r = size of each communities
% random sized k communities

% indexing each nodes of community
com_idx = cell(k,1);
prev = 0;
for i=1:k
    com_idx{i} = (1:n/k)+prev;
    prev = prev+n/k;
end

% using community label, generate adjacency matrix with parameter p and q
rm = tril(rand(n),-1)';
A = zeros(n);
for i=1:k
    for j=i:k
        sub_graph = rm(com_idx{i},com_idx{j});
        if i==j
            A(com_idx{i},com_idx{j}) = double(sub_graph>1-p);
        else
            A(com_idx{i},com_idx{j}) = double(sub_graph>1-q);
        end
    end
end

G = A+A';