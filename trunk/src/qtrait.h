/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate: 2006-12-23 23:02:07 -0600 (Sat, 23 Dec 2006) $
 *   $Rev: 644 $                                                      *
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

#ifndef _QTRAIT_H
#define _QTRAIT_H
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

	// ///////////////////////// Quantitative trait ////////////////////////

	/** \brief quantitative trait

	Genetic quantitative trait is tricky to simulate. In simuPOP, I employee
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

	class quanTrait: public Operator
	{
		public:
			/// constructor. default to be always active.
			quanTrait( int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr(1, "qtrait"))
				:Operator("","",stage, begin, end, step, at, rep, grp, infoFields)
			{
			}

			/// destructor
			virtual ~quanTrait()
			{
			}

			virtual Operator* clone() const
			{
				return new quanTrait(*this);
			}

			/// calculate/return quantitative trait etc
			virtual double qtrait(individual *)
			{
				///
				throw ValueError("This quantitative trait calculator is not supposed to be called directly");
				return 1.;
			}

			/// set qtrait to all individual
			bool apply(population& pop)
			{
				UINT idx = pop.infoIdx(infoField(0));
				// we need info to be in order
				GappedInfoIterator traitIt = pop.infoBegin(idx, true);
				for (population::IndIterator it = pop.indBegin(); it != pop.indEnd(); ++it)
					*traitIt++ = qtrait(&*it) ;

				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::qtrait::quantitative trait>" ;
			}
	};

	/** \brief quantitative trait according to genotype at one locus

	map selector. Assign qtrait value according to a
	given dictionary.
	*/

	class mapQuanTrait: public quanTrait
	{
		public:
			/** \brief create a map selector (quantitative trait according to genotype at one locus

			\param locus the locus index. The genotype of this locus will be axamed.
			\param loci the loci.
			\param qtrait a dictionary of qtrait. The genotype must be in the form of 'a-b'. This is the mean
			of quantitative trait. The actual trait value will be N(mean, sigma^2)
			For multiple loci, the form is 'a-b|c-d|e-f' etc.
			\param sigma standard deviation of the environmental facotr N(0,sigma^2).
			\param phase if true, a/b and b/a will have different qtrait value. Default to false.
			\param output and other parameters please refer to help(baseOperator.__init__)
			*/
			mapQuanTrait( vectoru loci, const strDict& qtrait, double sigma=0, bool phase=false,
				int stage=PostMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(1, "qtrait")):
			quanTrait(stage, begin, end, step, at, rep, grp, infoFields),
				m_loci(loci), m_dict(qtrait), m_sigma(sigma), m_phase(phase)
			{
			};

			virtual ~mapQuanTrait()
			{
			}

			virtual Operator* clone() const
			{
				return new mapQuanTrait(*this);
			}

			/// currently assuming diploid
			virtual double qtrait(individual * ind)
			{
				string key;

				for(vectoru::iterator loc=m_loci.begin(); loc!=m_loci.end(); ++loc)
				{
					/// get genotype of ind
					Allele a = ind->allele(*loc, 0);
					Allele b = ind->allele(*loc, 1);

					if( loc != m_loci.begin() )
						key += '|';

					if( ! m_phase && a > b )	  // ab=ba
						key +=  toStr(static_cast<int>(b)) + "-" + toStr(static_cast<int>(a));
					else
						key +=  toStr(static_cast<int>(a)) + "-" + toStr(static_cast<int>(b));
				}
				strDict::iterator pos = m_dict.find(key);

				DBG_ASSERT( pos != m_dict.end(), ValueError,
					"No qtrait value for genotype " + key);

				return( rng().randNormal(pos->second, m_sigma) );
			}

			virtual string __repr__()
			{
				return "<simuPOP::qtrait::map quantitative trait>" ;
			}

		private:
			/// one locus
			vectoru m_loci;

			/// qtrait for each genotype
			strDict m_dict;

			///
			double m_sigma;

			///
			bool m_phase;
	};

	/** \brief quantitative trait according to genotype at one locus

	multiple allele selector. This selector group alleles to disease
	and wild type and return qtrait to AA,Aa,aa. (A is wildtype).
	*/

	class maQuanTrait: public quanTrait
	{
		public:
			/** \brief create a multiple allele selector (quantitative trait according to diseased or wildtype
			alleles)

			\param locus the locus index. The genotype of this locus will be axamed.
			\param qtrait an array of qtrait of AA,Aa,aa. A is the wild type group.
			\param sigma an array of standard deviation for each of the trait genotype (AA, Aa, aa)
			\param wildtype an array of alleles in the wildtype group. Anything else is disease allele.
			   default = [0]
			\param output and other parameters please refer to help(baseOperator.__init__)
			*/
			maQuanTrait( vectoru loci, const vectorf& qtrait, const vectora& wildtype,
				const vectorf& sigma = vectorf(),
				int stage=PostMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(1, "qtrait")):
			quanTrait(stage, begin, end, step, at, rep, grp, infoFields),
				m_loci(loci), m_qtrait(qtrait), m_sigma(sigma), m_wildtype(wildtype)
			{
				if( m_sigma.empty())
					m_sigma.resize(3,0.);

				DBG_ASSERT( m_qtrait.size() == static_cast<UINT>(pow(static_cast<double>(3), 
					static_cast<double>(loci.size())))
					&& m_sigma.size() == m_qtrait.size(),
					ValueError, "Please specify qtrait for every combination of genotype.");
			};

			/// destructor
			virtual ~maQuanTrait()
			{
			}

			virtual Operator* clone() const
			{
				return new maQuanTrait(*this);
			}

			/// currently assuming diploid
			virtual double qtrait(individual * ind)
			{
				UINT index = 0;
				for(vectoru::iterator loc=m_loci.begin(); loc!=m_loci.end(); ++loc)
				{
					/// get genotype of ind
					Allele a = ind->allele(*loc, 0);
					Allele b = ind->allele(*loc, 1);

					int numWildtype=0;

					// count number of wildtype
					if( find(m_wildtype.begin(), m_wildtype.end(), a) != m_wildtype.end() )
						numWildtype ++;

					if( find(m_wildtype.begin(), m_wildtype.end(), b) != m_wildtype.end() )
						numWildtype ++;
					index = index*3 + 2-numWildtype;
				}

				return rng().randNormal(m_qtrait[index], m_sigma[index]);
			}

			virtual string __repr__()
			{
				return "<simuPOP::qtrait::multiple-alleles qtrait>" ;
			}

		private:
			/// one locus
			vectoru m_loci;

			/// qtrait for each genotype
			vectorf m_qtrait;

			///
			vectorf m_sigma;

			///
			vectora m_wildtype;
	};

	/** \brief quantitative trait according to genotype at multiple loci multiplicative model

	 multiple loci selector. This selector takes several selectors and
	 multiply their qtrait values...
	 e.g.
	   mlmquanTrait( [mapquanTrait(...), maquanTrait(...) ])
	 */

	class mlQuanTrait: public quanTrait
	{
		public:

#define QT_Multiplicative 1
#define QT_Additive 2

			/// vector of operator pointers.
			typedef std::vector< Operator * > vectorop;

		public:
			/** \brief multiple loci selector using a multiplicative model.

			\param qtraits a list of qtraits.
			*/
			mlQuanTrait( const vectorop qtraits, int mode = QT_Multiplicative, double sigma=0,
				int stage=PostMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(1, "qtrait")):
			quanTrait(stage, begin, end, step, at, rep, grp, infoFields),
				m_qtraits(0), m_sigma(sigma), m_mode(mode)
			{
				DBG_FAILIF( qtraits.empty(), ValueError, "Please specify at least one selector.");
				for(vectorop::const_iterator s = qtraits.begin(), sEnd=qtraits.end(); s != sEnd; ++s)
				{
					DBG_ASSERT( (*s)->__repr__().substr(10,6) == "qtrait", ValueError,
						"Expecting a vector of quantitative trait calculator");
					m_qtraits.push_back( (*s)->clone() );
				}

			};

			virtual ~mlQuanTrait()
			{
				for(vectorop::iterator s = m_qtraits.begin(), sEnd=m_qtraits.end(); s != sEnd; ++s)
					delete *s;
			}

			virtual Operator* clone() const
			{
				throw ValueError("Multi-loci selector can not be nested.");
			}

			/// currently assuming diploid
			virtual double qtrait(individual * ind)
			{
				if(m_mode == QT_Multiplicative )
				{
					double fit = 1;
					for(vectorop::iterator s = m_qtraits.begin(), sEnd=m_qtraits.end();
						s != sEnd; ++s)
					fit *= static_cast<quanTrait *>(*s)->qtrait(ind);
					return rng().randNormal(fit, m_sigma);
				}
				else if(m_mode == QT_Additive )
				{
					double fit = 0;
					for(vectorop::iterator s = m_qtraits.begin(), sEnd=m_qtraits.end();
						s != sEnd; ++s)
					fit +=  static_cast<quanTrait *>(*s)->qtrait(ind);
					return rng().randNormal(fit, m_sigma);
				}
				return 0.;
			}

			virtual string __repr__()
			{
				return "<simuPOP::qtrait::multiple-loci qtrait>" ;
			}

		private:
			/// a list of qtraits
			vectorop m_qtraits;

			///
			double m_sigma;

			/// mode
			int m_mode;
	};

	/** \brief quantitative trait using user supplied function

	Assign qtrait value by calling a user supplied function
	*/

	class pyQuanTrait: public quanTrait
	{
		public:
			/** \brief create a python hybrid selector

			\param loci susceptibility loci. The genotype at these loci will be
			passed to func.
			\param func a Python function that accept genotypes at susceptibility loci
			and return qtrait value.
			\param output and other parameters please refer to help(baseOperator.__init__)
			*/
			/// provide locus and qtrait for 11, 12, 13 (in the form of dictionary)
			pyQuanTrait( vectoru loci, PyObject* func,
				int stage=PostMating, int begin=0, int end=-1, int step=1,
				vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(1, "qtrait")):
			quanTrait(stage, begin, end, step, at, rep, grp, infoFields),
				m_loci(loci), m_alleles(0), m_len(0), m_numArray(NULL)
			{
				if( !PyCallable_Check(func))
					throw ValueError("Passed variable is not a callable python function.");

				Py_XINCREF(func);
				m_func = func;

				DBG_FAILIF( loci.empty(), ValueError,
					"Please specify susceptibility loci");
			};

			virtual ~pyQuanTrait()
			{
				if( m_func != NULL)
					Py_DECREF(m_func);
			}

			/// CPPONLY
			pyQuanTrait(const pyQuanTrait& rhs):
			quanTrait(rhs),
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
				return new pyQuanTrait(*this);
			}

			/// currently assuming diploid
			virtual double qtrait(individual * ind)
			{
				if( m_len == 0)
				{
					m_len = m_loci.size() * ind->ploidy();
					m_alleles.resize( m_len);
					m_numArray = Allele_Vec_As_NumArray( m_alleles.begin(), m_alleles.end() );
				}

				DBG_FAILIF( static_cast<size_t>(m_len) != ind->ploidy() * m_loci.size(),
					SystemError,
					"Length of m_len is wrong. Have you changed pop type?" );

				UINT pEnd = ind->ploidy();
				for(size_t i=0, iEnd=m_loci.size(), j=0; i < iEnd; ++i)
					for(UINT p=0; p < pEnd; ++p)
						m_alleles[j++] = ind->allele(m_loci[i], p);

				double resDouble;
				PyCallFunc(m_func, "(O)", m_numArray, resDouble, PyObj_As_Double);
				return resDouble;
			}

			virtual string __repr__()
			{
				return "<simuPOP::qtrait::python qtrait>" ;
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
