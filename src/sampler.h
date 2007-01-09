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

#ifndef _SAMPLER_H
#define _SAMPLER_H
/**
\file
\brief head file of class selector:public Operator
*/
#include "utility.h"
#include "operator.h"

#include "boost/tuple/tuple.hpp"
#include <numeric>
using std::min;

const string ASC_AS_Fields[2] = {"father_idx", "mother_idx"};

namespace simuPOP
{
	// ///////////////////// SUBSET ///////////////////////////////////////////
	// // ascertainment ........................................

	/// thrink population accroding to some outside value

	class pySubset: public Operator
	{

		public:
			/// create a directmigrator
			/**
			\param keep a carray of the length of population.
			its values will be assigned to info.
			\stage and other parameters please see help(baseOperator.__init__)
			*/
			pySubset(const vectori& keep=vectori(),
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()) :
			Operator( "", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_keep(keep)
			{
			}

			/// destructor
			virtual ~pySubset()
			{
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new pySubset(*this);
			}

			virtual bool apply(population& pop)
			{
				DBG_ASSERT( m_keep.size() >= pop.popSize() ,
					ValueError, "Given subpopid array has a length of "
					+ toStr( m_keep.size()) + " which is less than population size "
					+ toStr(pop.popSize()));

				for(size_t i=0, iEnd=pop.popSize(); i<iEnd; ++i)
				{
					DBG_ASSERT(static_cast<size_t>(m_keep[i]) <= MaxSubPopID, ValueError,
						"Subpop id exceeding maximum allowed subpopulations");
					// subpop id is short
					pop.ind(i).setSubPopID(static_cast<SubPopID>(m_keep[i]));
				}

				pop.setSubPopByIndID();
				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::pySubset>" ;
			}

		private:
			vectori m_keep;
	};

	/// sample from population and save samples

	/// sample operator will generate a new subpopulation in pop namespace.

	class sample: public Operator
	{

		public:
			/// create a sample
			/**
			\param name variable name of the sampled population (will be put in pop local namespace)
			\param nameExpr expression version of name. If both name and nameExpr is empty, do not store pop.
			\param times how many times to run the sample process? This is usually one, but we
			   may want to take several random samples.
			\param saveAs filename to save the population.
			\param saveAsExpr expression for save filename
			\param format to save sample(s)
			\param stage and other parameters please see help(baseOperator.__init__)
			*/
			sample( const string& name="sample", const string& nameExpr="", UINT times=1,
				const string& saveAs="", const string& saveAsExpr="",   const string& format="auto",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				: Operator( "", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_name(name), m_nameExpr(nameExpr,""), m_times(times), m_saveAs(saveAs),
				m_saveAsExpr(saveAsExpr), m_format(format)
			{
				if(name=="" && nameExpr=="" && saveAs=="" && saveAsExpr=="")
					throw ValueError("Please specify name (or nameExpr) or saveAs (or saveAsExpr) to save sample.");
			}

			/// destructor
			virtual ~sample(){};

			/// this function is very important
			virtual Operator* clone() const
			{
				return new sample(*this);
			}

			virtual bool prepareSample(population& )
			{
				return true;
			}

			virtual population& drawsample(population& pop)
			{
				// keep pop untouched.
				// use ind.info() directly
				throw SystemError("This function is not supposed to be called directly");
				return pop;
			}

			/// return the samples
			PyObject* samples(population& pop);

			virtual bool apply(population& pop);

			virtual string __repr__()
			{
				return "<simuPOP::sample>" ;
			}

		public:
			/// service functions that can be used by everyone

			/// save the idx of each individual to a filed
			/// usually 'oldindex'
			void saveIndIndex(population& pop, const string& indexField="oldindex");

			/// reset father_idx and mother_idx
			/// using the samed
			void resetParentalIndex(population& pop, const string& fatherField="father_idx",
				const string& motherField="mother_idx", const string& indexField="oldindex");

			/// find offspring and spouse
			void findOffspringAndSpouse(population& pop, unsigned ancestralDepth, unsigned maxOffspring,
				const string& fatherField, const string& motherField,
				const string& spouseField, const string& offspringField);

			/// set all sub pop id to -1 (remove)
			void resetSubPopID(population& pop);
		private:

			/// name to save sample, default to 'sample'
			string m_name;

			/// pop name
			Expression m_nameExpr;

			/// sample times
			UINT m_times;

			/// filename to save sample
			string  m_saveAs;

			/// saveas expression
			Expression m_saveAsExpr;

			/// format to save samples
			string m_format;
	};

	/// thrink population accroding to some outside value

	class randomSample: public sample
	{

		public:
			/// draw random sample, regardless of affected status
			/**
			\param size size of sample. It can be either a number, representing
			  the overall sample size, regardless of population strucutre;
			  or an array, representing number of samples drawn from each subpopulation.
			\param stage and other parameters please see help(baseOperator.__init__)
			\param name variable name of the sampled population (will be put in pop local namespace)
			\param nameExpr expression version of name. If both name and nameExpr is empty, do not store pop.
			\param times how many times to run the sample process? This is usually one, but we
			   may want to take several random samples.
			\param saveAs filename to save the population.
			\param saveAsExpr expression for save filename
			\param format to save sample(s)
			\param stage and other parameters please see help(baseOperator.__init__)

			\note ancestral populations will not be copied to the samples
			*/
			randomSample( vectorlu size=vectorlu(),
				const string& name="sample", const string& nameExpr="", UINT times=1,
				const string& saveAs="", const string& saveAsExpr="",   const string& format="auto",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				: sample(name, nameExpr, times, saveAs, saveAsExpr, format,
				stage, begin, end, step, at, rep, grp, infoFields),
				m_size(size)
			{
			}

			/// destructor
			virtual ~randomSample(){};

			/// this function is very important
			virtual Operator* clone() const
			{
				return new randomSample(*this);
			}

			/// value checking
			virtual bool prepareSample(population& pop);

			virtual population& drawsample(population& pop);

			virtual string __repr__()
			{
				return "<simuPOP::random sample>" ;
			}

		private:
			/// sample size
			vectorlu m_size;
	};

	/// thrink population accroding to some outside value

	class caseControlSample: public sample
	{

		public:
			/// draw cases and controls
			/**
			\param cases number of cases, or an array of number of cases from
			  each subpopulation.
			\param controls number of controls, or an array of number of controls
			  from each subpopulation.
			\param name variable name of the sampled population (will be put in pop local namespace)
			\param nameExpr expression version of name. If both name and nameExpr is empty, do not store pop.
			\param times how many times to run the sample process? This is usually one, but we
			   may want to take several random samples.
			\param saveAs filename to save the population.
			\param saveAsExpr expression for save filename
			\param format to save sample(s)
			\param stage and other parameters please see help(baseOperator.__init__)
			*/
			caseControlSample( const vectori& cases=vectori(), const vectori& controls = vectori(),
				bool spSample=false, const string& name="sample", const string& nameExpr="", UINT times=1,
				const string& saveAs="", const string& saveAsExpr="",   const string& format="auto",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				: sample(name, nameExpr, times, saveAs, saveAsExpr, format,
				stage, begin, end, step, at, rep, grp, infoFields),
				m_numCases(cases), m_numControls(controls), m_spSample(spSample),
				m_caseIdx(0), m_controlIdx(0)
			{
			}

			/// destructor
			virtual ~caseControlSample(){};

			/// this function is very important
			virtual Operator* clone() const
			{
				return new caseControlSample(*this);
			}

			virtual bool prepareSample(population& pop );

			virtual population& drawsample(population& pop);

			virtual string __repr__()
			{
				return "<simuPOP::case control sample>" ;
			}

		private:

			/// number of cases, use vectori instead of vectorlu because
			/// this will be post to setIntVectorVar
			vectori m_numCases, m_numControls;

			/// whether or not sample from each subpop
			bool m_spSample;

			vector< vectori > m_caseIdx, m_controlIdx;
	};

	/// thrink population accroding to some outside value

	class affectedSibpairSample: public sample
	{

		public:
			/// draw cases and controls
			/**
			\param size number of affected sibpairs to be sampled.
			  Can be a number or an array. If a number is given, it is
			  the total number of sibpairs, ignoring population structure.
			  Otherwise, given number of sibpairs are sampled from
			  subpopulations. If size is unspecified, this operator
			  will return all affected sibpairs.
			\param countOnly set variables about number of affected sibpairs,
			do not actually draw the sample
			\param name variable name of the sampled population (will be put in pop local namespace)
			\param nameExpr expression version of name. If both name and nameExpr is empty, do not store pop.
			\param times how many times to run the sample process? This is usually one, but we
			may want to take several random samples.
			\param saveAs filename to save the population.
			\param saveAsExpr expression for save filename
			\param format to save sample(s)
			\param stage and other parameters please see help(baseOperator.__init__)
			*/
			affectedSibpairSample(vectoru size = vectoru(),
				bool chooseUnaffected=false,
				bool countOnly=false,
				const string& name="sample", const string& nameExpr="", UINT times=1,
				const string& saveAs="", const string& saveAsExpr="",
				const string& format="auto",
				int stage=PostMating, int begin=0, int end=-1,
				int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(ASC_AS_Fields, ASC_AS_Fields+2))
				: sample(name, nameExpr, times, saveAs, saveAsExpr, format,
				stage, begin, end, step, at, rep, grp, infoFields),
				m_size(size), m_affectedness(!chooseUnaffected), m_countOnly(countOnly),
				m_validSibs(0)
			{
			}

			/// destructor
			virtual ~affectedSibpairSample(){};

			/// this function is very important
			virtual Operator* clone() const
			{
				return new affectedSibpairSample(*this);
			}

			virtual bool prepareSample(population& pop);

			virtual population& drawsample(population& pop);

			virtual string __repr__()
			{
				return "<simuPOP::affected sibpair sample>" ;
			}

		private:
			/// sample size
			vectoru m_size;

			bool m_affectedness;

			// do not draw sample
			bool m_countOnly;

			/// sibs for all subpopulations
			vector<vectorlu> m_validSibs;

			/// index to fields
			UINT m_father_id, m_mother_id;

	};

	class largePedigreeSample: public sample
	{

		public:
			/// draw cases and controls
			/**
			\param size number of affected sibpairs to be sampled.
			  Can be a number or an array. If a number is given, it is
			  the total number of sibpairs, ignoring population structure.
			  Otherwise, given number of sibpairs are sampled from
			  subpopulations. If size is unspecified, this operator
			  will return all affected sibpairs.
			\param countOnly set variables about number of affected sibpairs,
			do not actually draw the sample
			\param name variable name of the sampled population (will be put in pop local namespace)
			\param nameExpr expression version of name. If both name and nameExpr is empty, do not store pop.
			\param times how many times to run the sample process? This is usually one, but we
			may want to take several random samples.
			\param saveAs filename to save the population.
			\param saveAsExpr expression for save filename
			\param format to save sample(s)
			\param stage and other parameters please see help(baseOperator.__init__)
			*/
			largePedigreeSample(vectoru size = vectoru(),
				unsigned minTotalSize=0,
				unsigned maxOffspring=5,
				unsigned minPedSize=5, 
                unsigned minAffected=0,
				bool countOnly=false,
				const string& name="sample", const string& nameExpr="", UINT times=1,
				const string& saveAs="", const string& saveAsExpr="",
				const string& format="auto",
				int stage=PostMating, int begin=0, int end=-1,
				int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL,
				const vectorstr& infoFields=vectorstr(ASC_AS_Fields, ASC_AS_Fields+2))
				: sample(name, nameExpr, times, saveAs, saveAsExpr, format,
				stage, begin, end, step, at, rep, grp, infoFields),
				m_size(size), m_minTotalSize(minTotalSize), m_maxOffspring(maxOffspring),
				m_minPedSize(minPedSize), m_minAffected(minAffected),
				m_countOnly(countOnly), m_validPedigrees()
			{
			}

			/// destructor
			virtual ~largePedigreeSample(){};

			/// this function is very important
			virtual Operator* clone() const
			{
				return new largePedigreeSample(*this);
			}

			virtual bool prepareSample(population& pop);

			virtual population& drawsample(population& pop);

			virtual string __repr__()
			{
				return "<simuPOP::affected sibpair sample>" ;
			}

		private:
			typedef vector<boost::tuple<double, int> > pedArray;

			/// sample size
			vectoru m_size;

			/// control total size
			unsigned m_minTotalSize;

			///
			unsigned m_maxOffspring;

			///
			unsigned m_minPedSize;

			///
			unsigned m_minAffected;

			// do not draw sample
			bool m_countOnly;

			/// sibs for all subpopulations
			/// we need to also save size information.
			vector<pedArray> m_validPedigrees;
	};

	/// thrink population accroding to some outside value

	class pySample: public sample
	{

		public:
			/// create a python sampler
			/**
			\param keep a carray of the length of population.
			its values will be assigned to info.
			\param keepAncestralPop: -1 (all), 0 (no), 1(one ancestral pop) and so on.
			\stage and other parameters please see help(baseOperator.__init__)
			*/
			pySample( PyObject * keep, int keepAncestralPops=-1,
				const string& name="sample", const string& nameExpr="", UINT times=1,
				const string& saveAs="", const string& saveAsExpr="",   const string& format="auto",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr())
				: sample(name, nameExpr, times, saveAs, saveAsExpr, format,
				stage, begin, end, step, at, rep, grp, infoFields),
				m_keepAncestralPops(keepAncestralPops)
			{
				DBG_ASSERT( PyObj_Is_IntNumArray(keep), ValueError,
					"Passed vector is not a Python/Numeric int array");
				Py_INCREF(keep);
				m_keep = keep;
			}

			/// destructor
			virtual ~pySample()
			{
				if( m_keep != NULL)
					Py_DECREF(m_keep);
			}

			/// CPPONLY
			pySample(const pySample& rhs):
			sample(rhs),
				m_keep(rhs.m_keep)
			{
				if( m_keep != NULL)
					Py_INCREF(m_keep);
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new pySample(*this);
			}

			virtual population& drawsample(population& pop)
			{
				DBG_ASSERT( NumArray_Size(m_keep) >= static_cast<int>(pop.popSize()) ,
					ValueError, "Given subpopid array has a length of "
					+ toStr( NumArray_Size(m_keep)) + " which is less than population size "
					+ toStr(pop.popSize()));

				long * id = reinterpret_cast<long*>(NumArray_Data(m_keep));

				for(size_t i=0, iEnd=pop.popSize(); i<iEnd; ++i)
				{
					DBG_ASSERT(static_cast<size_t>(id[i]) <= MaxSubPopID, ValueError,
						"Subpop id exceeding maximum allowed subpopulations");
					// convert from int to signed short
					pop.ind(i).setSubPopID(static_cast<SubPopID>(id[i]));
				}

				return pop.newPopByIndID(m_keepAncestralPops);
			}

			virtual string __repr__()
			{
				return "<simuPOP::pySubset>" ;
			}

		private:
			PyObject* m_keep;

			int m_keepAncestralPops;
	};

}
#endif
