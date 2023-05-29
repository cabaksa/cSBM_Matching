function [AM] = M2adj(M)

m = size(M,1);
AM = zeros(m);
for i=1:m
    AM(M(i,1),M(i,2))=1;
end
