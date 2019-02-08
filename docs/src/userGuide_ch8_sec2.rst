Demographic model
=================

The original paper used a very simple instant population growth model. Under the
model assumption, a population with an initial population size :math:`N_{0}`
would evolve :math:`G_{0}` generations, instantly expand its population size to
:math:`N_{1}` and evolve another :math:`G_{1}` generations. Such a model can be
easily implemented as follows:

::

   def ins_expansion(gen):
       'An instant population growth model'
       if gen < G0:
           return N0
       else:
           return N1

Other demographic models could be implemented similarly. For example, an
exponential population growth model that expand the population size from
:math:`N_{0}` to :math:`N_{1}` in :math:`G_{1}` generations could be defined as

::

   def exp_expansion(gen):
       'An exponential population growth model'
       if gen < G0:
           return N0
       else:
           rate = (math.log(N1) - math.log(N0))/G1
           return int(N0 * math.exp((gen - G0) * rate))

That is to say, we first solve :math:`r` from
:math:`N_{1}=N_{0}\exp\left(rG_{1}\right)` and then calculate
:math:`N_{t}=N_{0}\exp\left(rG\right)` for a given generation.

There is a problem here: the above definitions treat ``N0``, ``G0``, ``N1`` and
``G1`` as global variables. This is OK for small scripts but is certainly not a
good idea for larger scripts especially when different parameters will be used.
A better way is to wrap these functions by another function that accept ``N0``,
``G0``, ``N1`` and ``G1`` as parameters. That is demonstrated in Example
:ref:`reichDemo <reichDemo>` where a function ``demo_model`` is defined to
return either an instant or an exponential population growth demographic
function.

.. _reichDemo:

**Example**: *A demographic function producer*

::

   >>> import simuPOP as sim
   >>> import math
   >>> def demo_model(model, N0=1000, N1=100000, G0=500, G1=500):
   ...     '''Return a demographic function 
   ...     model: linear or exponential
   ...     N0:   Initial sim.population size.
   ...     N1:   Ending sim.population size.
   ...     G0:   Length of burn-in stage.
   ...     G1:   Length of sim.population expansion stage.
   ...     '''
   ...     def ins_expansion(gen):
   ...         if gen < G0:
   ...             return N0
   ...         else:
   ...             return N1
   ...     rate = (math.log(N1) - math.log(N0))/G1
   ...     def exp_expansion(gen):
   ...         if gen < G0:
   ...             return N0
   ...         else:            
   ...             return int(N0 * math.exp((gen - G0) * rate))
   ...     if model == 'instant':
   ...         return ins_expansion
   ...     elif model == 'exponential':
   ...         return exp_expansion
   ... 
   >>> # when needed, create a demographic function as follows
   >>> demo_func = demo_model('exponential', 1000, 100000, 500, 500)
   >>> # sim.population size at generation 700
   >>> print(demo_func(700))
   6309

   now exiting runScriptInteractively...

`Download reichDemo.py <reichDemo.py>`_

.. note::

   The defined demographic functions return the total population size (a number) at
   each generation beacuse no subpopulation is considered. A list of subpopulation
   sizes should be returned if there are more than one subpopulations.


