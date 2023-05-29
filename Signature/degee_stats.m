function [comA_fdeg,comB_fdeg,avg_degA,var_degA,avg_degB,var_degB] = degee_stats(k,comA_deg,comB_deg)

comA_fdeg = cell(k,k);
comB_fdeg = cell(k,k);
avg_degA = zeros(k);
avg_degB = zeros(k);
var_degA = zeros(k);
var_degB = zeros(k);
for i=1:k
    for j=1:k
        avg_degA(i,j) = median(comA_deg{i,j});
        avg_degB(i,j) = median(comB_deg{i,j});
        var_degA(i,j) = var(comA_deg{i,j});
        var_degB(i,j) = var(comB_deg{i,j});
    end
end


for i=1:k
    for j=1:k
        comA_fdeg{i,j} = comA_deg{i,j}>avg_degA(i,j);
        comB_fdeg{i,j} = comB_deg{i,j}>avg_degB(i,j);
    end
end
