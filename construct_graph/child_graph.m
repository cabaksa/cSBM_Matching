function child = child_graph(pgraph,s)

m = size(pgraph,1);
int_sample = tril(rand(m)<s,-1);
sample = int_sample+int_sample.';
child = pgraph.*sample;