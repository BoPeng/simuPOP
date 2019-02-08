Utility Classes
===============


class WithArgs
--------------

.. class:: WithArgs

   This class wraps around a user-provided function and provides an
   attribute ``args`` so that simuPOP knows which parameters to send to the
   function. This is only needed if the function can not be defined with
   allowed parameters.

   .. method:: WithArgs.WithArgs(func, args)

      Return a callable object that wraps around function ``func``.
      Parameter ``args`` should be a list of parameter names.



class WithMode
--------------

.. class:: WithMode

   This class wraps around a user-provided output string, function
   or file handle (acceptable by parameter ``output`` of operators) so
   that simuPOP knows which mode the output should be written to. For
   example, if the output of the operator is a binary compressed stream,
   ``WithMode(output, 'b')`` could be used to tell the operators to 
   output bytes instead of string. This is most needed for Python 3 
   because files in Python 2 accepts string even if they are opened in
   binary mode.

   .. method:: WithMode.WithMode(output, mode='')

      Return an object that wraps around ``output`` and tells simuPOP
      to output string in ``mode``. This class currently only support
      ``mode=''`` for text mode and ``mode='b'`` for binary output.



class RNG
---------

.. class:: RNG

   This random number generator class wraps around a number of random
   number generators from GNU Scientific Library. You can obtain and
   change the  RNG used by the current  simuPOP module through the
   ``getRNG()`` function, or create a separate random number generator
   and use it in your script.


   .. method:: RNG(name=None, seed=0)


      Create a  RNG object using specified name and seed. If *rng* is
      not given, environmental variable ``GSL_RNG_TYPE`` will be used
      if it is available. Otherwise, generator ``mt19937`` will be
      used. If *seed* is not given, ``/dev/urandom``, ``/dev/random``,
      or other system random number source will be used to guarantee
      that random seeds are used even if more than one  simuPOP
      sessions are started simultaneously. Names of supported random
      number generators are available from
      ``moduleInfo()['availableRNGs']``.


   .. method:: RNG.name()

      Return the name of the current random number generator.

   .. method:: RNG.randBinomial(n, p)

      Generate a random number following a binomial distribution with
      parameters *n* and *p*.

   .. method:: RNG.randChisq(nu)

      Generate a random number following a Chi-squared distribution
      with *nu* degrees of freedom.

   .. method:: RNG.randExponential(mu)

      Generate a random number following a exponential distribution
      with parameter *mu*.

   .. method:: RNG.randGamma(a, b)

      Generate a random number following a gamma distribution with a
      shape parameters *a* and scale parameter *b*.

   .. method:: RNG.randGeometric(p)

      Generate a random number following a geometric distribution with
      parameter *p*.

   .. method:: RNG.randInt(n)

      return a random number in the range of ``[0, 1, 2, ... n-1]``

   .. method:: RNG.randMultinomial(N, p)

      Generate a random number following a multinomial distribution
      with parameters *N* and *p* (a list of probabilities).

   .. method:: RNG.randNormal(mu, sigma)

      Generate a random number following a normal distribution with
      mean *mu* and standard deviation *sigma*.

   .. method:: RNG.randPoisson(mu)

      Generate a random number following a Poisson distribution with
      parameter *mu*.

   .. method:: RNG.randTruncatedBinomial(n, p)

      Generate a positive random number following a zero-truncated
      binomial distribution with parameters *n* and *p*.

   .. method:: RNG.randTruncatedPoisson(mu)

      Generate a positive random number following a zero-truncated
      Poisson distribution with parameter *mu*.

   .. method:: RNG.randUniform()

      Generate a random number following a rng_uniform [0, 1)
      distribution.

   .. method:: RNG.seed()

      Return the seed used to initialize the  RNG. This can be used to
      repeat a previous session.

   .. method:: RNG.set(name=None, seed=0)

      Replace the existing random number generator using  RNG*name*
      with seed *seed*. If *seed* is 0, a random seed will be used. If
      *name* is empty, use the existing  RNG but reset the seed.


class WeightedSampler
---------------------

.. class:: WeightedSampler

   A random number generator that returns ``0``, ``1``, ..., ``k-1``
   with probabilites that are proportional to their weights. For
   example, a weighted sampler with weights ``4``, ``3``, ``2`` and
   ``1`` will return numbers ``0``, ``1``, ``2`` and ``3`` with
   probabilities ``0.4``, ``0.3``, ``0.2`` and ``0.1``, respectively.
   If an additional parameter ``N`` is specified, the weighted sampler
   will return exact proportions of numbers if ``N`` numbers are
   returned. The version without additional parameter is similar to
   the ``sample(prob, replace=FALSE)`` function of the R statistical
   package.


   .. method:: WeightedSampler(weights=[], N=0)


      Creates a weighted sampler that returns ``0``, ``1``, ...
      ``k-1`` when a list of ``k`` weights are specified (*weights*).
      *weights* do not have to add up to 1. If a non-zero *N* is
      specified, exact proportions of numbers will be returned in *N*
      returned numbers.


   .. method:: WeightedSampler.draw()

      Returns a random number between ``0`` and ``k-1`` with
      probabilities that are proportional to specified weights.

   .. method:: WeightedSampler.drawSamples(n=1)

      Returns a list of *n* random numbers


