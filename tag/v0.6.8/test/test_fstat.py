# testing fstat read/write/Fst calculation

from simuUtil import *
simu = simulator(population(subPop=[500,100,1000], ploidy=2,
   loci = [1], maxAllele=2), randomMating(), rep=2)

simu.apply( [
  initByValue(
    value = [[1,1],[1,2],[2,2],[1,1],[1,2],[2,2],[1,2],[1,2],[2,2]],
    indRange = [[0,125], [125,375],[375,500],[500,550],
    [550,580],[580,600],[600,700],[700, 1200], [1200,1600]])
 ] )

# now, write to a file
SaveFstat(simu.population(0), "p0.dat")
print open("p0.dat").read()
pop = LoadFstat("p0.dat")
dumper().apply(pop)

#
count = stat(
  popSize = 1,
  numOfMale = 1,
  alleleFreq = [0],
  heteroFreq = [0],
  genoFreq = [0],
  Fst = [ 0]
  )

simu.apply([count])
listVars(simu.dvars(0))
