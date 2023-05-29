function [A,B,com_idx] = generate_graph(n,p,q,k,s,perm)

[G,com_idx] = balanced_sbm(n,p,q,k);
%[G,com_idx] = unbalanced_sbm(n,p,q,k);
pG = permutation_matrix(G,n,perm);

A = child_graph(G,s);
B = child_graph(pG,s);