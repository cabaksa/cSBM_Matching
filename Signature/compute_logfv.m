function [fA,vA] = compute_logfv(To,comA_deg,c,ck,l)

%To = signature_tree(A,com_idxA,comA_fdeg,c,k1);

m = size(comA_deg{c,c},1);
log_com = log(comA_deg{c,c});
log_com(comA_deg{c,c}==0)=0;
Ec = median(log_com);
Vc = var(log_com);
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
        fA{i}(vec) = sum(log_com(To{i}{vec}))-sT*(Ec+1);
    end
end




