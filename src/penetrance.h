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
\brief head file of class selector:public Operator
*/
#include "utility.h"
#include "operator.h"

#include "boost/tuple/tuple.hpp"
#include <numeric>
using std::min;

namespace simuPOP
{
	// ///////////////////////// PENETRANCE ///////////////////////////////

	/** \brief penetrance

	Please refer to the user's guide for details.
	*/

	class penetrance: public Operator
	{
		public:
			/// constructor. default to be always active.
			// ancestralGen: 0: apply to current gen
			//               -1: apply to all gen
			//               otherwise, apply to n ancestral generations.
			/// If one field is specified, it will be used to store penetrance
			/// values.
			/// default to post mating
			penetrance(int ancestralGen=-1, int stage=DuringMating,
				int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr())
				:Operator("","",stage, begin, end, step, at, rep, grp, infoFields),
				m_ancestralGen(ancestralGen)
			{
			}

			/// destructor
			virtual ~penetrance()
			{
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new penetrance(*this);
			}

			/// calculate/return penetrance etc
			virtual double penet(individual *)
			{
				throw ValueError("This penetrance calculator is not supposed to be called directly");
				return 1.;
			}

			/// set pentrance to all individuals and record penetrance if requested.
			virtual bool apply(population& pop);

			/// set penetrance to all individual
			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL);

			virtual string __repr__()
			{
				return "<simuPOP::penetrance>" ;
			}

		private:
			/// how to handle ancestral gen
			int m_ancestralGen;
	};

	/** \brief penetrance according to genotype at one locus

	map selector. Assign penetrance value according to a
	given dictionary.
	*/

	class mapPenetrance: public penetrance
	{
		public:
			/** \brief create a map penetrance function (penetrance according to genotype at one locus

			\param locus the locus index. The genotype of this locus will be axamed.
			\param loci the loci index. The genotype of this locus will be axamed.
			\param penetrance a dictionary of penetrance. The genotype must be in the form of 'a-b' for single
			   locus.
			\param phase if true, a/b and b/a will have different penetrance value. Default to false.
			\param output and other parameters please refer to help(baseOperator.__init__)
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

			virtual Operator* clone() const
			{
				return new mapPenetrance(*this);
			}

			/// currently assuming diploid
			virtual double penet(individual * ind);

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

	/** \brief penetrance according to genotype at one locus

	multiple allele selector. This selector group alleles to disease
	and wild type and return penetrance to AA,Aa,aa. (A is wildtype).
	*/

	class maPenetrance: public penetrance
	{
		public:
			/** \brief create a multiple allele selector (penetrance according to diseased or wildtype
			alleles)

			\param locus the locus index. The genotype of this locus will be axamed.
			\param loci the loci index.
			\param penetrance an array of penetrance of AA,Aa,aa. A is the
			   wild type group. In the case of multiple loci, fitness should be in the order of
					BB Bb bb
				 AA 1  2  3
				 Aa 4  5  6
				 aa 7  8  9
			\param wildtype an array of alleles in the wildtype group. Anything else is disease allele.,
			default = [0]
			\param output and other parameters please refer to help(baseOperator.__init__)
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

			virtual Operator* clone() const
			{
				return new maPenetrance(*this);
			}

			/// currently assuming diploid
			virtual double penet(individual * ind);

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

	/** \brief penetrance according to genotype at multiple loci multiplicative model

	 multiple loci selector. This selector takes several selectors and
	 multiply their penetrance values...
	 e.g.
	   mlmpenetrance( [mappenetrance(...), mapenetrance(...) ])
	 */

	class mlPenetrance: public penetrance
	{
		public:

#define PEN_Multiplicative 1
#define PEN_Additive 2
#define PEN_Heterogeneity 3

			typedef std::vector< Operator * > vectorop;

		public:
			/** \brief multiple loci selector using a multiplicative model.

			\param selectors a list of selectors.
			\param mode one of PEN_Multiplicative, PEN_Additive, PEN_Heterogeneity
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

			virtual Operator* clone() const
			{
				throw ValueError("Multi-loci selector can not be nested.");
			}

			/// currently assuming diploid
			virtual double penet(individual * ind);

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

	/** \brief penetrance using user supplied function

	Assign penetrance value by calling a user supplied function
	*/

	class pyPenetrance: public penetrance
	{
		public:
			/** \brief create a python hybrid selector

			\param loci susceptibility loci. The genotype at these loci will be
			passed to func.
			\param func a Python function that accept genotypes at susceptibility loci
			and return penetrance value.
			\param output and other parameters please refer to help(baseOperator.__init__)
			*/
			/// provide locus and penetrance for 11, 12, 13 (in the form of dictionary)
			pyPenetrance( vectoru loci, PyObject* func, int ancestralGen=-1,
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

			virtual Operator* clone() const
			{
				return new pyPenetrance(*this);
			}

			/// currently assuming diploid
			virtual double penet(individual * ind);

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
