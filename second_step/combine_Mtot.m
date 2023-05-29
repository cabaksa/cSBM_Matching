function [Mtot] = combine_Mtot(perm_mat,com,pcom)

k = size(perm_mat,2);
M = cell(1,k);

for c=1:k
    m_mat = perm_mat{c};
    m = size(m_mat,1);
    M{c} = zeros(m,2);
    for i=1:m
        if sum(m_mat(i,:))==0
            continue
        end
        M{c}(i,1) = i;
        M{c}(i,2) = find(m_mat(i,:)==1);
    end
    M{c}(~any(M{c},2),:)=[];
end

Mtot =[];
for c=1:k
    x = com{c}(M{c}(:,1))';
    y = pcom{c}(M{c}(:,2))';
    new = [x,y];
    Mtot = [Mtot;new];
end