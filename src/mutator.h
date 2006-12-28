/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$                                                     *
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

#ifndef _MUTATOR_H
#define _MUTATOR_H
/**
\file
\brief head file of class mutator:public Operator
*/
/// for hybrid mutator
#include "operator.h"

namespace simuPOP
{
	/** \brief mutator class.

	Do not use this class directly. It just provide interface for real mutators.

	Every mutator can specify rate (equal rate) or rates (different rate for different
	loci) and a vector of applicable loci (default to all but should have the same
	length with rates if rates have length greater than one).

	max allele can be specified as well but more parameter, if needed, should
	be implemented by individual mutator classes.

	Number of possible allelic states: Most theoretical studies assume an infinite
	number of allelic states to avoid any homoplasy. If it facilitates analysis,
	this is however extremely unrealistic.

	@author Bo Peng
	*/

	class mutator: public Operator
	{
		public:
			/** \brief create a mutator
			All mutators have the following common parameters. However, the actual meaning
			of these parameters may vary according to different model. Check the manual
			for details!!! (help(kamMutator) for example.)

			\param rate single rate for all applicable loci (atLoci). Will be ignored
			if rates is specified; or it can be an array of rates, the same length as atLoci.
			\param atLoci a vector of loci index. Can be ignored only when single
			rate is specified. Default to all loci.
			\param maxAllele max allowable allele. Interpreted by each sub mutaor class. Default to pop.maxAllele().
			*/
			mutator( vectorf rate=vectorf(),
				vectori atLoci=vectori(),
				UINT maxAllele=0,
				string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				:Operator(output, outputExpr, stage, begin, end, step, at, rep, grp, infoFields),
				m_rate(rate), m_maxAllele(maxAllele), m_atLoci(atLoci),
				m_bt(rng()), m_initialized(false), m_mutCount(0)
			{
				if( m_rate.empty() )
					throw ValueError("You should specify a rate, or a sequence of rate.");

				if( rate.size() > 1 && atLoci.empty())
					throw ValueError("If you use variable rates, you should specify atLoci for each of the rate.");

				if( rate.size() > 1 && !atLoci.empty() && rate.size() != atLoci.size() )
					throw ValueError("If both rates and atLoci are specified, they should have the same length.");

#ifdef BINARYALLELE
				DBG_WARNING( maxAllele > 1, "MaxAllele for binary libraries must be 1");
				m_maxAllele = 1;
#else
				DBG_ASSERT( maxAllele <= MaxAllele, ValueError,
					"The maximum allele number exceeds " + toStr(MaxAllele)
					+ ". \nIf you need longer allele size, please use simuPOP_la libraries.");
#endif
			}

			/// destructor
			virtual ~mutator()
			{
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new mutator(*this);
			}

			/// return mutation rate
			vectorf rate()
			{
				return m_rate;
			}

			/// set an array of rates
			void setRate(const vectorf rate, const vectori atLoci = vectori())
			{
				if( rate.size() != 1 && rate.size() != atLoci.size() )
					throw ValueError("If you specify more than one rate values, you should also specify corresponding applicable loci");

				m_rate = rate;
				if( ! atLoci.empty())
					m_atLoci = atLoci;

				m_initialized = false;
			}

			/// return max allowable allele number
			UINT maxAllele()
			{
				return m_maxAllele;
			}

			///
			void setMaxAllele(UINT maxAllele)
			{
#ifndef BINARYALLELE
				m_maxAllele = maxAllele;
#endif
			}

			/// return mutation count
			ULONG mutationCount(size_t locus)
			{
				DBG_ASSERT( locus < m_mutCount.size(), IndexError,
					"locus index " + toStr(locus) + " is out of range");
				return m_mutCount[locus];
			}

			/// return mutation counts
			vectoru mutationCounts()
			{
				return m_mutCount;
			}

			/// how to mutate a single allele.
			/// this is usually the only function that need to be defined by the subclasses.
			virtual void mutate(AlleleRef allele)
			{
				throw SystemError("You are not supposed to call this base mutator funciton.");
			};

			/// apply!
			virtual bool apply(population& pop);

		private:

			/// initialize bernulli trial according to pop size etc
			virtual void initialize(population& pop);

		private:
			/// mutation rates
			vectorf m_rate;

			/// maxAllele
			UINT m_maxAllele;

			/// applicable loci.
			vectori m_atLoci;

			/// bernulli trials. bitSet mutation results.
			BernulliTrials m_bt;

			/// initialized? the first apply() call will trigger an initialization process.
			bool m_initialized;

			/// report the number of mutation events
			vectoru m_mutCount;
	};

	/// K-Allele Model mutator
	/**
	Under this model, there are K (here refers as maxAllele) possible allele states, and any
	allele has a constant probability (rate/(K-1)) of mutation towards any of the K-1 allelic states.

	\note the theoretical mutation rate is rates/(K-1)  towards any of the K-1 allelic states.
	So rates is actually the probability to mutate!

	\sa Crow & Kimura 1970
	*/

	class kamMutator: public mutator
	{
		public:
			/// K-Allele Model mutator
			/**
			\param rate  mutation rate. It is 'probability to mutate'. The actual
			   mutation rate to any of the other K-1 allelic states are rates/(K-1)!
			\param atLoci and other parameters: refer to help(mutator), help(baseOperator.__init__)
			\param maxAllele maxAllele that can be mutated to. For binary libraries
			  allelic states will be [0, maxAllele]. For others, they are [1, maxAllele]
			*/
			kamMutator(const vectorf& rate=vectorf(),
				const vectori& atLoci=vectori(),
				UINT maxAllele=0,
				string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				:mutator( rate, atLoci, maxAllele,
				output, outputExpr, stage, begin, end, step, at, rep, grp, infoFields)
			{
			}

			~kamMutator(){}

			/// mutate to a state other than current state with equal probability
			virtual void mutate(AlleleRef allele);

			/// this function is very important
			virtual Operator* clone() const
			{
				return new kamMutator(*this);
			}

			virtual string __repr__()
			{
				return "<simuPOP::k-allele model mutator K=" +
					toStr(this->maxAllele()) + ">" ;
			}
	};

	/// stepwise mutation model.
	/**
	Stepwise mutation model (SMM) assumes that alleles are represented by integer values
	and that a mutation either increases or decreases the allele value by one.
	For variable number tandem repeats loci (VNTR), the allele value is generally
	taken as the number of tandem repeats in the DNA sequence.

	\sa Kimura & Ohta 1978
	*/

	class smmMutator: public mutator
	{
		public:
			///
			/**
			The stepwise mutation model (SMM) is developed for allozymes. It  provides better description
			 for these kinds of evolutionary processes.

			\param rate: mutation rate
			\param incProb probability to increase allele state. Default to 1
			\param atLoci and other parameters: refer to help(mutator), help(baseOperator.__init__)

			*/
			smmMutator(vectorf rate=vectorf(), vectori atLoci=vectori(),
				UINT maxAllele=0, double incProb=0.5,
				string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				:mutator( rate, atLoci, maxAllele,
				output, outputExpr, stage, begin, end, step, at, rep, grp, infoFields),
				m_incProb(incProb)
			{
#ifdef BINARYALLELE
				DBG_WARNING(true, "Symetric stepwise mutation does not work well on two state alleles.");
#endif
				DBG_ASSERT( fcmp_ge( incProb, 0.) && fcmp_le( incProb, 1.),
					ValueError, "Inc probability should be between [0,1], given " + toStr(incProb));
			}

			~smmMutator(){}

			virtual void mutate(AlleleRef allele)
			{
				// inc
				if( rng().randUniform01() < m_incProb)
				{
					if( AlleleUnsigned(allele) < this->maxAllele() )
						AlleleInc(allele);
				}
				// dec (use !=0 instead of > 0 to avoid warning inbinary mode
				else if( allele != 0 )
					AlleleDec(allele);
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new smmMutator(*this);
			}

			virtual string __repr__()
			{
				return "<simuPOP::step-wise mutation model mutator>" ;
			}

		private:
			/// probability to increase allele state
			double m_incProb;
	};

	/// stepwise mutation model.
	/**
	Generalized Stepwise mutation model (GSM) assumes that alleles are represented by integer values
	and that a mutation either increases or decreases the allele value by a random value.

	\sa Kimura & Ohta 1978
	*/

	class gsmMutator: public mutator
	{
		public:
			///
			/**
			The generalized stepwise mutation model (GMM) is developed for allozymes.
			It  provides better description for these kinds of evolutionary processes.

			\param rate: mutation rate
			\param incProb probability to increase allele state. Default to 0.5
			\param atLoci and other parameters: refer to help(mutator), help(baseOperator.__init__)
			\param func return number of steps. no parameter
			*/
			gsmMutator( vectorf rate=vectorf(), vectori atLoci=vectori(),
				UINT maxAllele=0, double incProb=0.5, double p=0, PyObject* func=NULL,
				string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				:mutator( rate, atLoci, maxAllele,
				output, outputExpr, stage, begin, end, step, at, rep, grp, infoFields),
				m_incProb(incProb), m_p(p), m_func(func)
			{
				DBG_ASSERT( fcmp_ge( incProb, 0.) && fcmp_le( incProb, 1.),
					ValueError, "Inc probability should be between [0,1], given " + toStr(incProb));

#ifdef BINARYALLELE
				DBG_WARNING(true, "Generalized stepwise mutation does not work well on two state alleles.");
#endif

				if( func != NULL)				  // use this function
				{
					DBG_ASSERT( PyCallable_Check(func),
						ValueError, "Func is not a Python function");

					Py_XINCREF(func);
					m_func = func;
				}
				else
				{
					DBG_ASSERT( fcmp_ge( p, 0.) && fcmp_le( p, 1.),
						ValueError, "Parameter p of a geometric distribution should be between [0,1], given " + toStr(m_p));
				}
			}

			~gsmMutator(){}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new gsmMutator(*this);
			}

			virtual void mutate(AlleleRef allele);

			virtual string __repr__()
			{
				return "<simuPOP::generalized step-wise mutator>" ;
			}

		private:
			/// probability to increase allele state
			double m_incProb;

			/// parameter for geometric gsm
			double m_p;

			/// the function to return random number
			PyObject* m_func;
	};

	/// mixed mutation model . has not been implemented.

	class pyMutator: public mutator
	{
		public:
			pyMutator(vectorf rate=vectorf(), vectori atLoci=vectori(), UINT maxAllele=0,
				PyObject* func=NULL,
				string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				:mutator( rate,atLoci, maxAllele,
				output, outputExpr, stage, begin, end, step, at, rep, grp, infoFields), m_func(NULL)
			{
				if( !PyCallable_Check(func))
					throw ValueError("Passed variable is not a callable python function.");

				Py_XINCREF(func);
				m_func = func;
			}

			~pyMutator()
			{
				if( m_func != NULL )
					Py_DECREF(m_func);
			}

			/// CPPONLY
			pyMutator(const pyMutator& rhs):
			mutator(rhs), m_func(rhs.m_func)
			{
				if( m_func != NULL )
					Py_INCREF(m_func);
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new pyMutator(*this);
			}

			virtual void mutate(AlleleRef allele)
			{
				int resInt;
				PyCallFunc(m_func, "(i)", static_cast<unsigned int>(allele), resInt, PyObj_As_Int);
#ifdef BINARYALLELE                
                DBG_ASSERT(resInt == 0 || resInt == 1, ValueError,
                    "Can only mutate to 0 or 1 in binary mode.");
				allele = resInt != 0;
#else
                DBG_ASSERT(resInt <= MaxAllele, ValueError,
                    "Mutated to an allele greater than maximum allowed allele value");
				allele = static_cast<Allele>(resInt);
#endif                
			}

			virtual string __repr__()
			{
				return "<simuPOP::python mutator>" ;
			}

		private:

			PyObject* m_func;
	};

	/// \brief point mutator
	/** mutate specified individuals at specified loci to spcified allele.
	I.e., this is a non-random mutator used to introduce disease etc.
	*/

	class pointMutator: public Operator
	{
		public:
			/** \brief mutate once

			\param atLoci a vector of loci index.
			\param inds mutate 'inds' individuals
			\param toAllele mutate to 'toAllele'
			*/
			pointMutator(
				const vectori& atLoci,
				Allele toAllele,
				vectoru atPloidy=vectoru(),
				vectorlu inds=vectorlu(),
				string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				:Operator(output, outputExpr, stage, begin, end, step, at, rep, grp, infoFields),
				m_atLoci(atLoci), m_toAllele(toAllele),
				m_atPloidy(atPloidy), m_inds(inds), m_mutCount(0)
			{
				if(m_atPloidy.empty())
					m_atPloidy.push_back(0);
			}

			/// destructor
			virtual ~pointMutator()
			{
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new pointMutator(*this);
			}

			/// apply!
			virtual bool apply(population& pop)
			{
				m_mutCount.resize(pop.totNumLoci(), 0);
				// mutate each mutable locus
				for( size_t i=0, iEnd=m_atLoci.size(); i < iEnd; ++i)
				{
					for( vectorlu::iterator ind = m_inds.begin();
						ind != m_inds.end(); ++ind)
					{
						for( size_t p=0; p<m_atPloidy.size(); ++p)
						{
							m_mutCount[m_atLoci[i]]++;
							*(pop.ind(*ind).genoBegin(m_atPloidy[p])+m_atLoci[i])= m_toAllele;
						}
					}
				}								  // each applicable loci

				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::point mutator>" ;
			}

			/// return mutation count
			ULONG mutationCount(size_t locus)
			{
				DBG_ASSERT( locus < m_mutCount.size(), IndexError,
					"locus index " + toStr(locus) + " is out of range");
				return m_mutCount[locus];
			}

			/// return mutation counts
			vectoru mutationCounts()
			{
				return m_mutCount;
			}

		private:

			/// applicable loci.
			vectori m_atLoci;
			Allele m_toAllele;
			vectoru m_atPloidy;
			vectorlu m_inds;
			/// report the number of mutation events
			vectoru m_mutCount;
	};

}
#endif
