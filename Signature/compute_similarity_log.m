function [sim] = compute_similarity_log(A,B,com_idxA,com_idxB,c,l,sV,k,k1)

[comA_deg,comB_deg] = community_degree(A,B,com_idxA,com_idxB);
[comA_fdeg,comB_fdeg,~,~,~,~] = degee_stats(k,comA_deg,comB_deg);
TA = signature_tree(A,com_idxA,comA_fdeg,c,l,k1);
TB = signature_tree(B,com_idxB,comB_fdeg,c,l,k1);

% EcA = avg_degA(c,c);
% EcB = avg_degB(c,c);
% VcA = var_degA(c,c);
% VcB = var_degB(c,c);

m = size(com_idxA{c},2);
ck = size(k1,2);

[fA,vA] = compute_logfv(TA,comA_deg,c,ck,l);
[fB,vB] = compute_logfv(TB,comB_deg,c,ck,l);

cV = sV;
sim = zeros(m);
for i=1:m
    for j=1:m
        x = (fA{i}(cV)-fB{j}(cV)).^2;
        y = vA{i}(cV)+vB{j}(cV);
        z = x./y;
        z(isnan(z)) = 0;
        sim(i,j) = sum(z);
    end
end

