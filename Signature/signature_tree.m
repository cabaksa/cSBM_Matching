function To = signature_tree(A,com_idx,com_fdeg,c,l,k1)

m = size(com_idx{c},2);
ga = graph(A(com_idx{c},com_idx{c}));
% neighbor infor
comA_nb = cell(m,1);
for i=1:m
    comA_nb{i} = neighbors(ga,i);
end

% compute l length sphere rooted by i.
sphere = cell(m,l);
for i=1:m
    dist = distances(ga,i);
    for d = 1:l
        sphere{i,d} = find(dist==d);
    end
end
% d = [];
% e = [];
% ek = [c,d,e];
% k1 = setdiff((1:k),ek);
[dgcode_set] = degree_code(com_fdeg,k1,c);
ck = size(k1,2);

To = cell(m,1);
for i=1:m
    To{i} = {i};
    for d = 0:(l-1)
        Tn = cell(2^(ck*(d+1)),1);
        for sd = 1:2^(ck*d)
            Tnb = unique(cell2mat(comA_nb(To{i}{sd})));
            lnode = intersect(sphere{i,d+1},Tnb);

            for su = 1:2^ck
                Tn{2^(ck)*(sd-1)+su,1} = intersect(lnode,dgcode_set{su});
            end
        end
        To{i} = Tn;
    end
end






