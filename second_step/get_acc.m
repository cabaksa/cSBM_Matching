function acc = get_acc(M,perm,com,pcom,c)

diff = perm(com{c}(M(:,1)))-pcom{c}(M(:,2));
acc = size(find(diff==0),2)/size(M,1);

