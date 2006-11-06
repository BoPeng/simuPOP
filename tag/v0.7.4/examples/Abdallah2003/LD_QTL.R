# deal with pvalues
pvalues = read.table('pvalues.log', colClasses='numeric', sep=',')
dim(pvalues)
p = apply(pvalues, 2 , mean)
plot(p, main='mean pvalues', xlab='marker distance', ylab='pvalue')
