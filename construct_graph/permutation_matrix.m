function pG = permutation_matrix(G,n,perm)

PM = full(sparse((1:n), perm, 1, n, n));
pG = PM'*G*PM;