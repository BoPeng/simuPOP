from simuPOP_ba import *
pop = population(10, loci=[2])
p = pop.arrGenotype()
print p
print p[2]
p[2] = 1
p[4] = 1
print p
Dump(pop)
