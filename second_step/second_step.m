function [Mtot] = second_step(app_adj,A,B,c,com,pcom,T)

Ac = A(com{c},com{c});
Bc = B(pcom{c},pcom{c});
[AM] = refinement_matching(Ac,Bc,app_adj,T);

k = size(com,1);
ek = setdiff((1:k),c);
Ek = size(ek,2);

perm_mat = cell(1,k);
perm_mat{c} = AM;

for i=1:Ek
    app_adj = seeded_matching(A,B,com,pcom,c,ek(i),AM);
    Ai = A(com{ek(i)},com{ek(i)});
    Bi = B(pcom{ek(i)},pcom{ek(i)});
    [perm_mat{ek(i)}] = refinement_matching(Ai,Bi,app_adj,T);
end
Mtot = combine_Mtot(perm_mat,com,pcom);
% acc1 = get_acc(Mtot);
