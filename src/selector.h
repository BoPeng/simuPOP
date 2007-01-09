/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$                                                      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef _SELECTOR_H
#define _SELECTOR_H
/**
\file
\brief head file of class selector:public Operator
*/
#include "utility.h"
#include "operator.h"

#include "boost/tuple/tuple.hpp"
#include <numeric>
using std::min;

namespace simuPOP
{
	/** \brief selection

	Genetic selection is tricky to simulate. In simuPOP, I employee
	an ability (fitness) to mate approach. Namely, the probability
	that an individual will be chosen for mating is proportional
	to its fitness value. More specifically,

	- PreMating selectors assign fitness values to each individual.

	- Sexless mating (e.g. binomialSelection) : individuals are chosen
	at probabilities that are proportional to their fitness values. More
	specifically, if there are N individuals with fitness values
	\f$f_i, i=1,...,N \f$, individual \f$i\f$ will have probability
	\f$ \frac{f_i}{\sum_{j=1}^N f_j} \f$ to be chosen to be passed
	to the next generation.

	- Random mating with sex (e.g. randommating): males and females are
	separated and each are chosen as described above.

	Please refer to the user's guide for details.
	*/

	class selector: public Operator
	{
		public:
			/// constructor. default to be always active.
			selector( int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr(1, "fitness"))
				:Operator("","",stage, begin, end, step, at, rep, grp, infoFields)
			{
			}

			/// destructor
			virtual ~selector()
			{
			}

			virtual Operator* clone() const
			{
				return new selector(*this);
			}

			/// calculate/return w11 etc
			virtual double indFitness(individual *)
			{
				///
				throw ValueError("This selector is not supposed to be called directly");
				return 1.;
			}

			/// set fitness to all individual
			bool apply(population& pop)
			{
				UINT fit_id = pop.infoIdx(this->infoField(0));
				GappedInfoIterator fitness = pop.infoBegin(fit_id, true);

				for (population::IndIterator it = pop.indBegin(); it != pop.indEnd(); ++it)
					*fitness++ = indFitness(&*it) ;

				// indicate selection is on.
				pop.setBoolVar("selection", true);
				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::selector>" ;
			}
	};

	/** \brief selection according to genotype at one locus

	map selector. Assign fitness value according to a
	given dictionary.
	*/

	class mapSelector: public selector
	{
		public:
			/** \brief create a map selector (selection according to genotype at one locus

			\param locus the locus index. The genotype of this locus will be axamed.
			\param loci the loci index. The genotype of this locus will be axamed.
			\param fitness a dictionary of fitness. The genotype must be in the form of 'a-b'
			   for single locus, and 'a-b|c-d|e-f' for multi-locus..
			\param phase if true, a/b and b/a will have different fitness value. Default to false.
			\param output and other parameters please refer to help(baseOperator.__init__)
			*/
			mapSelector( vectoru loci, const strDict& fitness, bool phase=false,
				int stage=PreMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(1, "fitness")):
			selector(stage, begin, end, step, at, rep, grp, infoFields),
				m_loci(loci), m_dict(fitness), m_phase(phase)
			{
			};

			virtual ~mapSelector()
			{
			}

			virtual Operator* clone() const
			{
				return new mapSelector(*this);
			}

			/// currently assuming diploid
			virtual double indFitness(individual * ind);

			virtual string __repr__()
			{
				return "<simuPOP::selector::map selector>" ;
			}

		private:
			/// one locus
			vectoru m_loci;

			/// fitness for each genotype
			strDict m_dict;

			///
			bool m_phase;
	};

	/** \brief selection according to genotype at one locus

	multiple allele selector. This selector group alleles to disease
	and wild type and return fitness to AA,Aa,aa. (A is wildtype).
	*/

	class maSelector: public selector
	{
		public:
			/** \brief create a multiple allele selector (selection according to diseased or wildtype
			alleles)
			Note that maSelector only work for diploid population now.

			\param locus the locus index. The genotype of this locus will be axamed.
			\param loci the loci index.
			\param fitness For the single locus case, fitness is an array of fitness of AA,Aa,aa. A is the
			   wild type group. In the case of multiple loci, fitness should be in the order of
				   BB Bb bb
				AA 1  2  3
				Aa 4  5  6
			aa 7  8  9
			The length for such table is 3^(#loci).
			\param wildtype an array of alleles in the wildtype group. Anything else is disease allele.
			default = [0]
			NOTE that wildtype at all loci are the same.
			\param output and other parameters please refer to help(baseOperator.__init__)
			*/
			maSelector( vectoru loci, const vectorf& fitness, const vectora& wildtype,
				int stage=PreMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(1, "fitness")):
			selector(stage, begin, end, step, at, rep, grp, infoFields),
				m_loci(loci), m_fitness(fitness), m_wildtype(wildtype)
			{
				DBG_ASSERT( m_fitness.size() == static_cast<UINT>(pow(static_cast<double>(3), 
					static_cast<double>(loci.size()))),
					ValueError, "Please specify fitness for each combination of genotype.");
			};

			virtual ~maSelector()
			{
			}

			virtual Operator* clone() const
			{
				return new maSelector(*this);
			}

			/// currently assuming diploid
			virtual double indFitness(individual * ind);

			virtual string __repr__()
			{
				return "<simuPOP::selector::multiple-alleles selector>" ;
			}

		private:
			/// one locus
			vectoru m_loci;

			/// fitness for each genotype
			vectorf m_fitness;

			///
			vectora m_wildtype;
	};

	/** \brief selection according to genotype at multiple loci multiplicative model

	 multiple loci selector. This selector takes several selectors and
	 multiply their fitness values...
	 e.g.
	   mlmselector( [mapselector(...), maselector(...) ])
	 */

	class mlSelector: public selector
	{
		public:

#define SEL_None 0
#define SEL_Multiplicative 1
#define SEL_Additive 2
#define SEL_Heterogeneity 3

			typedef std::vector< Operator * > vectorop;

		public:
			/** \brief multiple loci selector using a multiplicative model.

			\param selectors a list of selectors.
			*/
			mlSelector( const vectorop selectors, int mode = SEL_Multiplicative,
				int stage=PreMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(1, "fitness")):
			selector(stage, begin, end, step, at, rep, grp, infoFields),
				m_selectors(0), m_mode(mode)
			{
				DBG_FAILIF( selectors.empty(), ValueError, "Please specify at least one selector.");
				for(vectorop::const_iterator s = selectors.begin(), sEnd=selectors.end(); s != sEnd; ++s)
				{
					DBG_ASSERT( (*s)->__repr__().substr(10,8) == "selector", ValueError,
						"Expecting a list of fitness calculator. Given " + (*s)->__repr__());
					m_selectors.push_back( (*s)->clone());
				}
			};

			virtual ~mlSelector()
			{
				for(vectorop::iterator s = m_selectors.begin(), sEnd=m_selectors.end(); s != sEnd; ++s)
					delete *s;
			}

			virtual Operator* clone() const
			{
				throw ValueError("Multi-loci selector can not be nested.");
			}

			/// currently assuming diploid
			virtual double indFitness(individual * ind);

			virtual string __repr__()
			{
				return "<simuPOP::selector::multiple-loci selector>" ;
			}

		private:
			/// a list of selectors
			vectorop m_selectors;

			/// mode
			int m_mode;
	};

	/** \brief selection using user supplied function

	Assign fitness value by calling a user supplied function
	*/

	class pySelector: public selector
	{
		public:
			/** \brief create a python hybrid selector

			\param loci susceptibility loci. The genotype at these loci will be
			passed to func.
			\param func a Python function that accept genotypes at susceptibility loci
			and return fitness value.
			\param output and other parameters please refer to help(baseOperator.__init__)
			*/
			/// provide locus and fitness for 11, 12, 13 (in the form of dictionary)
			pySelector( vectoru loci, PyObject* func,
				int stage=PreMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(1, "fitness")):
			selector(stage, begin, end, step, at, rep, grp, infoFields),
				m_loci(loci), m_alleles(0), m_len(0), m_numArray(NULL)
			{
				if( !PyCallable_Check(func))
					throw ValueError("Passed variable is not a callable python function.");

				Py_XINCREF(func);
				m_func = func;

				DBG_FAILIF( loci.empty(), ValueError,
					"Please specify susceptibility loci");
			};

			/// destructor
			virtual ~pySelector()
			{
				if( m_func != NULL)
					Py_DECREF(m_func);
			}

			/// CPPONLY
			pySelector(const pySelector& rhs):
			selector(rhs),
				m_loci(rhs.m_loci),
				m_func(rhs.m_func),
				m_alleles(rhs.m_alleles),
				m_len(rhs.m_len),
				m_numArray(NULL)
			{
				if( m_func != NULL)
					Py_INCREF(m_func);
			}

			virtual Operator* clone() const
			{
				return new pySelector(*this);
			}

			/// currently assuming diploid
			virtual double indFitness(individual * ind);

			virtual string __repr__()
			{
				return "<simuPOP::selector::python selector>" ;
			}

		private:

			/// susceptibility loci
			vectoru m_loci;

			/// user supplied python function
			PyObject* m_func;

			/// copy of alleles of each individual a time.
			vectora m_alleles;

			/// length of m_alleles
			int m_len;

			/// the object that passed to func
			PyObject * m_numArray;

	};

}
#endif
