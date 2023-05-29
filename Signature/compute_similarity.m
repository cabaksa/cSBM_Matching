function [sim] = compute_similarity(A,B,com_idxA,com_idxB,c,l,sV,k,k1)

[comA_deg,comB_deg] = community_degree(A,B,com_idxA,com_idxB);
[comA_fdeg,comB_fdeg,avg_degA,var_degA,avg_degB,var_degB] = degee_stats(k,comA_deg,comB_deg);
TA = signature_tree(A,com_idxA,comA_fdeg,c,l,k1);
TB = signature_tree(B,com_idxB,comB_fdeg,c,l,k1);

EcA = avg_degA(c,c);
EcB = avg_degB(c,c);
VcA = var_degA(c,c);
VcB = var_degB(c,c);

m1 = size(com_idxA{c},2);
m2 = size(com_idxB{c},2);
ck = size(k1,2);

[fA,vA] = compute_fv(TA,comA_deg,c,ck,l,EcA,VcA);
[fB,vB] = compute_fv(TB,comB_deg,c,ck,l,EcB,VcB);

cV = sV;
sim = zeros(m1,m2);
for i=1:m1
    for j=1:m2
        x = (fA{i}(cV)-fB{j}(cV)).^2;
        y = vA{i}(cV)+vB{j}(cV);
        z = x./y;
        z(isnan(z)) = 0;
        sim(i,j) = sum(z);
    end
end