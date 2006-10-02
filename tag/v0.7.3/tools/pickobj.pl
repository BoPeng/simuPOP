#!/usr/bin/perl

# get the line

#

$cmdBefore= "g++ -shared -nostdlib  /usr/site/gcc-3.3/lib/gcc-lib/sparc-sun-solaris2.8/3.3/crti.o /usr/ccs/lib/values-Xa.o /usr/site/gcc-3.3/lib/gcc-lib/sparc-sun-solaris2.8/3.3/crtbegin.o  .libs/libsimuPOP_la-utility.o  .libs/libsimuPOP_la-borosh13.o .libs/libsimuPOP_la-fishman2x.o .libs/libsimuPOP_la-mt.o .libs/libsimuPOP_la-rand.o .libs/libsimuPOP_la-ranmar.o .libs/libsimuPOP_la-types.o .libs/libsimuPOP_la-cmrg.o .libs/libsimuPOP_la-gfsr4.o .libs/libsimuPOP_la-r250.o .libs/libsimuPOP_la-random.o .libs/libsimuPOP_la-rng.o .libs/libsimuPOP_la-uni32.o .libs/libsimuPOP_la-coveyou.o .libs/libsimuPOP_la-knuthran2.o .libs/libsimuPOP_la-ran0.o .libs/libsimuPOP_la-randu.o .libs/libsimuPOP_la-slatec.o .libs/libsimuPOP_la-uni.o .libs/libsimuPOP_la-default.o .libs/libsimuPOP_la-knuthran.o .libs/libsimuPOP_la-ran1.o .libs/libsimuPOP_la-ranf.o .libs/libsimuPOP_la-taus113.o .libs/libsimuPOP_la-vax.o .libs/libsimuPOP_la-file.o .libs/libsimuPOP_la-lecuyer21.o .libs/libsimuPOP_la-ran2.o .libs/libsimuPOP_la-ranlux.o .libs/libsimuPOP_la-taus.o .libs/libsimuPOP_la-waterman14.o .libs/libsimuPOP_la-fishman18.o .libs/libsimuPOP_la-minstd.o .libs/libsimuPOP_la-ran3.o .libs/libsimuPOP_la-ranlxd.o .libs/libsimuPOP_la-transputer.o .libs/libsimuPOP_la-zuf.o .libs/libsimuPOP_la-fishman20.o .libs/libsimuPOP_la-mrg.o .libs/libsimuPOP_la-rand48.o .libs/libsimuPOP_la-ranlxs.o .libs/libsimuPOP_la-tt.o";

$cmdAfter= " .libs/libsimuPOP_la-simuPOP_wrap.o  -Wl,-R -Wl,/usr/site/gcc-3.3/lib/. -Wl,-R -Wl,/usr/site/gcc-3.3/lib/. -L/build/gcc-3.3/src/obj/sunos5/gcc -L/build/gcc-3.3/src/obj/sunos5/sparc-sun-solaris2.8/libstdc++-v3/src/.libs -L/build/gcc-3.3/src/obj/sunos5/sparc-sun-solaris2.8/libstdc++-v3/src -L/usr/site/matlab/share/extern/lib/sol2 -lmx -lmat -leng -L/usr/site/gcc-3.3/lib/gcc-lib/sparc-sun-solaris2.8/3.3 -L/usr/ccs/bin -L/usr/ccs/lib -L/usr/site/gcc-3.3/lib/gcc-lib/sparc-sun-solaris2.8/3.3/../../.. /usr/site/gcc-3.3/lib/./libstdc++.so -lm -lgcc /usr/site/gcc-3.3/lib/gcc-lib/sparc-sun-solaris2.8/3.3/crtend.o /usr/site/gcc-3.3/lib/gcc-lib/sparc-sun-solaris2.8/3.3/crtn.o  -Wl,-h -Wl,libsimuPOP.so.0 -o .libs/libsimuPOP.so.0.0.0; cp -f .libs/libsimuPOP.so.0.0.0 _simuPOP.so;python -c 'import _simuPOP'";

$objLine =".libs/libsimuPOP_la-minmax.o .libs/libsimuPOP_la-infnan.o .libs/libsimuPOP_la-coerce.o .libs/libsimuPOP_la-fdiv.o .libs/libsimuPOP_la-pow_int.o .libs/libsimuPOP_la-psi.o .libs/libsimuPOP_la-trig.o .libs/libsimuPOP_la-exp.o .libs/libsimuPOP_la-log.o .libs/libsimuPOP_la-erfc.o .libs/libsimuPOP_la-zeta.o .libs/libsimuPOP_la-elementary.o .libs/libsimuPOP_la-gamma.o .libs/libsimuPOP_la-bernoulli.o .libs/libsimuPOP_la-erlang.o .libs/libsimuPOP_la-gumbel.o .libs/libsimuPOP_la-nbinomial.o .libs/libsimuPOP_la-beta.o .libs/libsimuPOP_la-exponential.o .libs/libsimuPOP_la-hyperg.o .libs/libsimuPOP_la-pareto.o .libs/libsimuPOP_la-bigauss.o .libs/libsimuPOP_la-exppow.o .libs/libsimuPOP_la-landau.o .libs/libsimuPOP_la-pascal.o .libs/libsimuPOP_la-binomial.o .libs/libsimuPOP_la-fdist.o .libs/libsimuPOP_la-laplace.o .libs/libsimuPOP_la-poisson.o .libs/libsimuPOP_la-binomial_tpe.o .libs/libsimuPOP_la-flat.o .libs/libsimuPOP_la-levy.o .libs/libsimuPOP_la-rayleigh.o .libs/libsimuPOP_la-cauchy.o .libs/libsimuPOP_la-logarithmic.o .libs/libsimuPOP_la-shuffle.o .libs/libsimuPOP_la-rdgamma.o .libs/libsimuPOP_la-chisq.o .libs/libsimuPOP_la-gauss.o .libs/libsimuPOP_la-logistic.o .libs/libsimuPOP_la-sphere.o .libs/libsimuPOP_la-dirichlet.o .libs/libsimuPOP_la-gausstail.o .libs/libsimuPOP_la-lognormal.o .libs/libsimuPOP_la-tdist.o .libs/libsimuPOP_la-discrete.o .libs/libsimuPOP_la-geometric.o .libs/libsimuPOP_la-multinomial.o .libs/libsimuPOP_la-weibull.o .libs/libsimuPOP_la-error.o";

# get objs
@objs = split(/ /, $objLine);

@allObj = @objs;

for ( $i=0 ; $i < @objs; $i++ ){
    # remove $objs[$i] from @allObj
    $curObj = $objs[$i];
    $cmd = $cmdBefore;
    @cmdObj = ();
    foreach $o (@allObj){
	
	if( $o ne $curObj){
	    push(@cmdObj, $o);
	    $cmd .= " " . $o . " ";
	}
    }
    $cmd .= $cmdAfter;
    print "dealing with $curObj \t";
  #  print "$cmd\n";
    if( system($cmd) == 0 ){ # success
      # remove curObj from @allObj
      print "Removing $curObj \n";
      @allObj = @cmdObj;
    }
  # exit;
}
foreach $a (@allObj){
  print "$a \n";
}



