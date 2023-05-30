clear
addpath("Signature")
addpath("construct_graph")
addpath("second_step")
addpath("baseline")

total_task = 20;
total_cor = 11;

% m : community size
% k : the number of community 
% p : edge connectivity of intra community
% q : edge connectivity of outra community
% c : comparing community
% l : length of tree
% sV : signature vector
% k1 : the number of communities which used in signature vector

m = 833;
k = 6;
n = m*k;
c = 1;
p = 0.025;
q = p/3;
perm = randperm(n);
GM = [(1:n)' (perm)'];
Ground_Truth = full(sparse(GM(:, 1), GM(:, 2), 1, n, n));
T = 5;
l = 2;
k1 = setdiff((1:k),[c 2]);
sV = (1:2^(size(k1,2)*l));

% Compute fraction of correctly matched pairs 
Fir_Our = zeros(total_task,total_cor);
Fir_GP1 = zeros(total_task,total_cor);
Fir_DP1 = zeros(total_task,total_cor);
Fir_GP2 = zeros(total_task,total_cor);
Fir_DP2 = zeros(total_task,total_cor);

Sec_Our = zeros(total_task,total_cor);
Sec_GP1 = zeros(total_task,total_cor);
Sec_DP1 = zeros(total_task,total_cor);
Sec_DP2 = zeros(total_task,total_cor);
Sec_GP2 = zeros(total_task,total_cor);

noise = zeros(1,total_cor);

%% re-code for community. 

for t = 1:total_task
    for alp = 1:total_cor
%% Construct correlated graphs 
        s = 1-(0.025*(alp-1));
        [A,B,com] = generate_graph(n,p,q,k,s,perm);
        pcom = cellfun(@(x) sort(perm(x)),com,'UniformOutput',false);
%% Our Algorithm
% <Input arguments>
% A,B : correlated adjacency matrix
% com,pcom = community label of nodes

        Sim_Our = compute_similarity(A,B,com,pcom,c,l,sV,k,k1);
        MP_Sig = matchpairs(-Sim_Our, -99999, 'max');
        Fir_Our(t,alp) = get_acc(MP_Sig,perm,com,pcom,c);
        
        Adj_MP_sig = M2adj(MP_Sig);
        P_sig = second_step(Adj_MP_sig,A,B,c,com,pcom,T);
        Sec_Our(t,alp) = size(find((P_sig(:,2)-perm(P_sig(:,1)))==0),1)/n;
%% DP 1
        Sim_ADP = matching_deg_pro(A,B);
        Sim_DP1 = Sim_ADP(com{c},pcom{c});
        MP_DP1 = matchpairs(Sim_DP1, -99999, 'max');
        Fir_DP1(t,alp) = get_acc(MP_DP1,perm,com,pcom,c);
        
        Adj_MDP1 = M2adj(MP_DP1);
        P_DP1 = second_step(Adj_MDP1,A,B,c,com,pcom,T);
        Sec_DP1(t,alp) = size(find((P_DP1(:,2)-perm(P_DP1(:,1)))==0),1)/n;
%% DP 2
        Fir_DP2_P = zeros(1,k);
        Acc_SecDP_c = zeros(1,k);
        for ci = 1:k
            Aci = A(com{ci},com{ci});
            Bci = B(pcom{ci},pcom{ci});
            Sim_DP2_c = matching_deg_pro(Aci,Bci);
            M_DP2_c = matchpairs(Sim_DP2_c, -99999, 'max');
            Fir_DP2_P(ci) = get_acc(M_DP2_c,perm,com,pcom,ci);
            
            Adj_MDP2_c = M2adj(M_DP2_c);
            Pc = refinement_matching(Aci,Bci,Adj_MDP2_c,T);
            Acc_SecDP_c(ci) = sum(diag(Pc*Ground_Truth(com{ci},pcom{ci})'));
        end
        
        Fir_DP2(t,alp) = mean(Fir_DP2_P);
        Sec_DP2(t,alp) = sum(Acc_SecDP_c)/n;
        
%% GP 1
        Sim_AGP = grampa(A,B,1);
        Sim_GP1 = Sim_AGP(com{c},pcom{c});        
        MP_GP1 = matchpairs(Sim_GP1, -99999, 'max');
        Adj_MGP1 = M2adj(MP_GP1);
        Fir_GP1(t,alp) = get_acc(MP_GP1,perm,com,pcom,c);
        P_GP1 = second_step(Adj_MGP1,A,B,c,com,pcom,T);
        Sec_GP1(t,alp) = size(find((P_GP1(:,2)-perm(P_GP1(:,1)))==0),1)/n;  
        
%% GP 2
        Fir_GP2_P = zeros(1,k);
        Acc_SecGP_c = zeros(1,k);
        
        for ci = 1:k
            Aci = A(com{ci},com{ci});
            Bci = B(pcom{ci},pcom{ci});
            Sim_GP2_c = grampa(Aci,Bci,1);
            M_GP2_c = matchpairs(Sim_GP2_c, -99999, 'max');
            Fir_GP2_P(ci) = get_acc(M_GP2_c,perm,com,pcom,ci);
            
            Adj_MGP2_c = M2adj(M_GP2_c);
            Pc = refinement_matching(Aci,Bci,Adj_MGP2_c,T);
            Acc_SecGP_c(ci) = sum(diag(Pc*Ground_Truth(com{ci},pcom{ci})'));
        end
        
        Fir_GP2(t,alp) = mean(Fir_GP2_P);
        Sec_GP2(t,alp) = sum(Acc_SecGP_c)/n;   

        noise(alp) = s;
    end
end
