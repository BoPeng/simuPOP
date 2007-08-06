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

#ifndef _PENETRANCE_H
#define _PENETRANCE_H
/**
\file
\brief head file of class selector:public baseOperator
*/
#include "utility.h"
#include "operator.h"

#include "boost/tuple/tuple.hpp"
#include <numeric>
using std::min;

namespace simuPOP
{
	// ///////////////////////// PENETRANCE ///////////////////////////////

	/// basic class of a penetrance operator
	/**
	Penetrance is the probability that one will have the disease when he has certain
	genotype(s). Calculation and the parameter set of penetrance are similar to those
	of fitness. An individual will be randomly marked as affected/unaffected according
	to his penetrance value.??? For example, an individual will have probability 0.8 to
	be affected if the penetrance is 0.8. \n

	Penetrance can be applied at any stage (default to \c DuringMating). It will be
	calculated during mating, and then the affected status will be set for each offspring.
	Penetrance can also be used as \c PreMating, \c PostMating or even \c PrePostMating??? operator.
	In these cases, the affected status will be set to all individuals according to their
	penetrance values. It is also possible to store penetrance in a given information
	field specified by \c infoFields parameter (e.g. <tt>infoFields=['penetrance']</tt>). This is
	useful to check the penetrance values at a later time. \n

	Affected status will be used for statistical purpose, and most importantly, ascertainment.
	They will be calculated along with fitness although they might not be used at every
	generation. You can use two operators: one for fitness/selection, active at every
	generation; one for affected status, active only at ascertainments, to avoid unnecessary
	calculation of the affected status. \n

	Pentrance values are used to set the affectedness of individuals, and are usually not saved.
	If you would like to know the penetrance value, you need to
	\li use <tt>addInfoField('penetrance')</tt> to the population to analyze. (Or use \c infoFields
	parameter of the population constructor), and
	\li use e.g., <tt>mlPenetrance(...., infoFields=['penetrance'])</tt> to add the penetrance field
	to the penetrance operator you use. You may choose a name other than \c 'penetrance' as long as
	the field names for the operator and population match.

	Penetrance functions can be applied to the current, all, or certain number of ancestral generations.
	This is controlled by the \c ancestralGen parameter, which is default to \c -1 (all available
	ancestral generations). You can set it to \c 0 if you only need affection??? status for the current
	generation, or specify a number \c n for the number of ancestral generations (n + 1 total generations)
	to process. Note that \c ancestralGen parameter is ignored if the penetrance operator is used
	as a during mating operator.
	*/
	class penetrance: public baseOperator
	{
		public:
			/// create a penetrance operator
			/**
			default to be always active.
			\param ancestralGen if this parameter is set to be \c 0, then apply penetrance to
				the current generation; if \c -1, apply to all generations; otherwise, apply
				to the specified number of ancestral generations
			\param stage specify the stage this operator will be applied, default to \c DuringMating.
			\param infoFields If one field is specified, it will be used to store penetrance values.???
			*/
			penetrance(int ancestralGen=-1, int stage=DuringMating,
				int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr())
				: baseOperator("","",stage, begin, end, step, at, rep, grp, infoFields),
				m_ancestralGen(ancestralGen)
			{
			}

			/// destructor
			virtual ~penetrance()
			{
			}

			/// deep copy of a penetrance operator
			virtual baseOperator * clone() const
			{
				return new penetrance(*this);
			}

			/// calculate/return penetrance etc.
			virtual double penet(individual *)
			{
				throw ValueError("This penetrance calculator is not supposed to be called directly");
				return 1.;
			}

			/// set penetrance to all individuals and record penetrance if requested
			virtual bool apply(population& pop);

			/// set penetrance to all individuals
			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL);

			/// used by Python print function to print out the general information of the penetrance operator
			virtual string __repr__()
			{
				return "<simuPOP::penetrance>" ;
			}

		private:
			/// how to handle ancestral gen
			int m_ancestralGen;
	};

	/// penetrance according to the genotype at one locus
	/**
	Assign penetrance using a table with keys \c 'X-Y' where \c X and \c Y are allele numbers.
	<funcForm>MapPenetrance</funcForm>
	*/
	class mapPenetrance: public penetrance
	{
		public:
			/// create a map penetrance operator
			/**
			\param locus the locus index. The genotype of this locus will be examed.???
			\param loci the loci indices. The genotypes of these loci will be examed.
			\param penetrance a dictionary of penetrance. The genotype must be in the form
				of 'a-b' for a single locus.
			\param phase if True, <tt>a/b</tt> and <tt>b/a</tt> will have different penetrance values.
				Default to \c False.
			\param output and other parameters please refer to help(baseOperator.__init__)???
			*/
			mapPenetrance( vectoru loci, const strDict& penet, bool phase=false,
				int ancestralGen=-1, int stage=DuringMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr()):
			penetrance(ancestralGen, stage, begin, end, step, at, rep, grp, infoFields),
				m_loci(loci), m_dict(penet), m_phase(phase)
			{
			};

			virtual ~mapPenetrance()
			{
			}

			/// deep copy of a map penetrance operator
			virtual baseOperator * clone() const
			{
				return new mapPenetrance(*this);
			}

			/// currently assuming diploid???
			virtual double penet(individual * ind);

			/// used by Python print function to print out the general information of the map penetrance operator
			virtual string __repr__()
			{
				return "<simuPOP::penetrance::map penetrance>" ;
			}

		private:
			/// one locus
			vectoru m_loci;

			/// penetrance for each genotype
			strDict m_dict;

			///
			bool m_phase;
	};

	/// multiple allele penetrance operator
	/**
	This is called 'multiple-alleles'??? penetrance. It separates alleles into two groups:
	wildtype and disease alleles. Wildtype alleles are specified by parameter \c wildtype
	and any other alleles are considered as diseased alleles. \c maPenetrance accepts an
	array of fitness for AA, Aa, aa in the single-locus case, and a longer table for
	multi-locus case. Penetrance is then set for any given genotype.
	<funcForm>MaPenetrance</funcForm>
	*/
	class maPenetrance: public penetrance
	{
		public:
			/// create a multiple allele penetrance operator (penetrance according to diseased or wildtype alleles)
			/**
			\param locus the locus index. The genotype of this locus will be examed.???
			\param loci the loci indices. The genotypes of these loci will be examed.
			\param penetrance an array of penetrance values of AA, Aa, aa. A is the
				wild type group. In the case of multiple loci, fitness should be in the order of
				AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb.
			\param wildtype an array of alleles in the wildtype group. Any other alleles will
				be considered as in the disease allele group.
			\param output and other parameters please refer to help(baseOperator.__init__)???
			*/
			maPenetrance( vectoru loci, const vectorf& penet, const vectora& wildtype,
				int ancestralGen=-1,
				int stage=DuringMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr()):
			penetrance(ancestralGen, stage, begin, end, step, at, rep, grp, infoFields),
				m_loci(loci), m_penetrance(penet), m_wildtype(wildtype)
			{
				DBG_ASSERT( m_penetrance.size() ==  static_cast<UINT>(pow(static_cast<double>(3),
					static_cast<double>(loci.size()))),
					ValueError, "Please specify penetrance for each combination of genotype.");
			};

			virtual ~maPenetrance()
			{
			}

			/// deep copy of a multi-allele penetrance operator
			virtual baseOperator * clone() const
			{
				return new maPenetrance(*this);
			}

			/// currently assuming diploid???
			virtual double penet(individual * ind);

			/// used by Python print function to print out the general information of the multi-allele penetrance operator
			virtual string __repr__()
			{
				return "<simuPOP::penetrance::multiple-alleles penetrance>" ;
			}

		private:
			/// one locus
			vectoru m_loci;

			/// penetrance for each genotype
			vectorf m_penetrance;

			///
			vectora m_wildtype;
	};
	/// penetrance according to the genotype according to a multiple loci multiplicative model
	/**
	\c mlPentrance is the 'multiple-loci'??? penetrnace calculator. It accepts a list of
	penetrances and combine them according to the \c mode parameter, which takes one of the
	following values:

	\li \c PEN_Multiplicative: the penetrance is calculated as \f$ f=\prod f_{i} \f$.
	\li \c PEN_Additive: the penetrance is calculated as \f$ f=\min\left(1,\sum f_{i}\right) \f$.
		\f$ f \f$ will be set to \c 1 when \f$ f<0 \f$. In this case, \f$ s_{i} \f$??? are
		added, not \f$ f_{i} \f$ directly.
	\li \c PEN_Heterogeneity: the penetrance is calculated as \f$ f=1-\prod\left(1-f_{i}\right) \f$.

	Please refer to Neil Risch (1990) for detailed information about these models.

	<funcForm>MlPenetrance</funcForm>
	*/
	class mlPenetrance: public penetrance
	{
		public:

#define PEN_Multiplicative 1
#define PEN_Additive 2
#define PEN_Heterogeneity 3

			typedef std::vector< baseOperator * > vectorop;

		public:
			/// create a multiple loci penetrance operator using a multiplicative model
			/**
			\param peneOps a list of selectors???
			\param mode can be one of \c PEN_Multiplicative, \c PEN_Additive, and \c PEN_Heterogeneity
			*/
			mlPenetrance( const vectorop peneOps, int mode = PEN_Multiplicative,
				int ancestralGen=-1, int stage=DuringMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr()):
			penetrance(ancestralGen, stage, begin, end, step, at, rep, grp, infoFields),
				m_peneOps(0), m_mode(mode)
			{
				DBG_FAILIF( peneOps.empty(), ValueError, "Please specify at least one penetrance operator.");
				for(vectorop::const_iterator s = peneOps.begin(), sEnd=peneOps.end(); s != sEnd; ++s)
				{
					DBG_ASSERT( (*s)->__repr__().substr(10,10) == "penetrance", ValueError,
						"Expecting a list of penetrance calculator");

					m_peneOps.push_back( (*s)->clone() );
				}

			};

			virtual ~mlPenetrance()
			{
				for(vectorop::iterator s = m_peneOps.begin(), sEnd=m_peneOps.end(); s != sEnd; ++s)
					delete *s;
			}

			/// deep copy of a multi-loci penetrance operator
			virtual baseOperator * clone() const
			{
				throw ValueError("Multi-loci selector can not be nested.");
			}

			/// currently assuming diploid???
			virtual double penet(individual * ind);

			/// used by Python print function to print out the general information of the multiple-loci penetrance operator
			virtual string __repr__()
			{
				return "<simuPOP::penetrance::multiple-loci penetrance>" ;
			}

		private:
			/// a list of peneOps
			vectorop m_peneOps;

			/// mode
			int m_mode;
	};

	/// assign penetrance values by calling a user provided function
	/**
	For each individual, users provide a function to calculate penetrance. This
	method is very flexible but will be slower than previous operators since a
	function will be called for each individual.
	<funcForm>PyPenetrance</funcForm>
	*/
	class pyPenetrance: public penetrance
	{
		public:
			/**
			\param loci disease susceptibility loci. The genotypes at these loci will
				be passed to the provided Python function in the form of <tt>loc1_1, loc1_2, loc2_1, loc2_2, ...</tt>
				if the individuals are diploid.
			\param func a user-defined Python function that accepts an array of genotypes
				at susceptibility loci and return a penetrance value. The returned value
				should be between \c 0 and \c 1.
			\param output and other parameters please refer to help(baseOperator.__init__)???
			*/
			/// provide locus and penetrance for 11, 12, 13 (in the form of dictionary)
			pyPenetrance(const vectoru & loci, PyObject* func, int ancestralGen=-1,
				int stage=DuringMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr()):
			penetrance(ancestralGen, stage, begin, end, step, at, rep, grp, infoFields),
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
			virtual ~pyPenetrance()
			{
				if( m_func != NULL)
					Py_DECREF(m_func);
				if( m_numArray != NULL)
					Py_DECREF(m_numArray);
			}

			/// CPPONLY
			pyPenetrance(const pyPenetrance& rhs):
			penetrance(rhs),
				m_loci(rhs.m_loci),
				m_func(rhs.m_func),
				m_alleles(rhs.m_alleles),
				m_len(rhs.m_len),
				m_numArray(NULL)
			{
				if( m_func != NULL)
					Py_INCREF(m_func);
			}

			/// deep copy of a Python penetrance operator
			virtual baseOperator * clone() const
			{
				return new pyPenetrance(*this);
			}

			/// currently assuming diploid???
			virtual double penet(individual * ind);

			/// used by Python print function to print out the general information of the Python penetrance operator
			virtual string __repr__()
			{
				return "<simuPOP::penetrance::python penetrance>" ;
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
