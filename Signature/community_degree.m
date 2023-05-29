function [comA_deg,comB_deg] = community_degree(A,B,com_idxA,com_idxB)

k = size(com_idxA,1);
comA_deg = cell(k,k);
comB_deg = cell(k,k);

for i=1:k
    for j=1:k
        Ga = A(com_idxA{i},com_idxA{j});
        Gb = B(com_idxB{i},com_idxB{j});
        comA_deg{i,j} = sum(Ga,2);
        comB_deg{i,j} = sum(Gb,2);
    end
end