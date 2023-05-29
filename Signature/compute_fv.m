function [fA,vA] = compute_fv(To,comA_deg,c,ck,l,Ec,Vc)

%To = signature_tree(A,com_idxA,comA_fdeg,c,k1);

m = size(comA_deg{c,c},1);
fA = cell(m,1);
vA = cell(m,1);
% Ec = avg_degA(c,c);
% Vc = var_degA(c,c);

tv = 2^(ck*l);

for i=1:m
    fA{i} = zeros(tv,1);
    vA{i} = zeros(tv,1);
    for vec = 1:tv
        sT = size(To{i}{vec},1);
        vA{i}(vec) = Vc*sT;
        fA{i}(vec) = sum(comA_deg{c,c}(To{i}{vec}))-sT*(Ec+1);
    end
end




