function app_per = seeded_matching(A,B,com_idxA,com_idxB,c,d,pos_adj)

a = A(com_idxA{d},com_idxA{c});
b = B(com_idxB{c},com_idxB{d});
t = a*pos_adj*b;

m1 = size(com_idxA{d},2);
m2 = size(com_idxB{d},2);

M = zeros(min(m1,m2),2);

for i=1:min(m1,m2)
    if sum(t(i,:))==0
        continue
    else
        M(i,1) = i;
        [~,M(i,2)] = max(t(i,:));
        t(i,:)=0;
        t(:,M(i,2))=0;
    end
end
M(~any(M,2),:)=[];

app_per = zeros(m1,m2);
app_pair = size(M,1);

for i=1:app_pair
    app_per(M(i,1),M(i,2)) = 1;
end