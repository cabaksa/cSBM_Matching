function [AM] = refinement_matching(A,B,X,T)

% input : A,B, similairty matrix X, iteration num T
% T = 100;
n = size(A,1);
nS = A*X*B;
for i=1:T
    nM = matchpairs(nS, -99999, 'max');
    eperm = zeros(n);
    for j=1:n
        eperm(nM(j,1),nM(j,2))=1;
    end
    nS = A*eperm*B;
end
% acc = get_acc(nM);
AM = M2adj(nM);