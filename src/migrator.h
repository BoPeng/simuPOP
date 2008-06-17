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

#ifndef _MIGRATOR_H
#define _MIGRATOR_H
/**
 \file
 \brief head file of class migrator:public baseOperator
 */
#include "operator.h"
#include <list>
using std::list;

#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
using std::sort;
using std::random_shuffle;

#include "simuPOP_cfg.h"
#include <math.h>
// using std::sqrt;

#include "virtualSubPop.h"

namespace simuPOP {

typedef std::vector<vsp> vectorvsp;

/// migrate individuals from a (sub)population to another (sub)population
/**
   Migrator is the only way to mix genotypes of several subpopulations
   because mating is strictly within subpopulations in simuPOP. Migrators
   are quite flexible in simuPOP in the sense that
 \li migration can happen from and to a subset of subpopulations.
 \li migration can be done by probability, proportion or by counts. In
   	the case of probability, if the migration rate from subpopulation
 \c a to \c b is \c r, then everyone in subpopulation \c a will have this
   	probability to migrate to \c b. In the case of proportion, exactly
   	<tt>r*size_of_subPop_a</tt> individuals (chosen by random) will migrate
   	to subpopulation \c b. In the last case, a given number of individuals will
   migrate.
 \li new subpopulation can be generated through migration. You simply
   need to migrate to a subpopulation with a new subpopulation number.
 */
class migrator : public baseOperator
{

public:
#define MigrByProbability  1
#define MigrByProportion   2
#define MigrByCounts       3

public:
	/// create a migrator
	/**
	 \param rate migration rate, can be a proportion or counted number. Determined by
	   	parameter \c mode. \c rate should be an m by n matrix. If a number is given,
	   	the migration rate will be a \c m by \c n matrix of value \c r
	 \param mode one of \c MigrByProbability (default), \c MigrByProportion or \c MigrByCounts
	 \param fromSubPop an array of 'from' (virtual) subpopulations. Default to all. 
		If a single (virtual) subpopulation is specified, <tt>[]</tt> can be ignored.
		A virtual subpopulation should be as <tt>vsp(subPop, virtualSubPop)</tt>. For example,
		if you define a virtual subpopulation by sex, you can use
		<tt>fromSubPop=vsp(0, 0)</tt> to choose migrants only from the first virtual subpopulation
		of subpopulation 0.
	 \param toSubPop an array of 'to' subpopulations. Default to all
		subpopulations. If a single subpopulation is specified,
		<tt>[]</tt> can be ignored.
	 \param stage default to \c PreMating

	 \note
	 \li The overall population size will not be changed. (Mating schemes can
	   do that). If you would like to keep the subpopulation sizes after migration, you
	   can use the \c newSubPopSize or \c newSubPopSizeExpr parameter of a mating scheme.
	 \li \c rate is a matrix with dimensions determined by \c fromSubPop and \c toSubPop.
	   By default, \c rate is a matrix with element \c r(i,j), where \c r(i, j) is the
	   migration rate, probability or count from subpopulation \c i to \c j. If \c fromSubPop
	   and/or \c toSubPop are given, migration will only happen between these subpopulations.
	   An extreme case is 'point migration', <tt>rate=[[r]], fromSubPop=a, toSubPop=b</tt>
	   which migrate from subpopulation \c a to \c b with given rate \c r.
	 */
	migrator(const matrix & rate, int mode = MigrByProbability,
	         const vectorvsp & fromSubPop = vectorvsp(), vectoru toSubPop = vectoru(),
	         int stage = PreMating, int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
	         int rep = REP_ALL, int grp = GRP_ALL, const vectorstr & infoFields = vectorstr())
		: baseOperator("", "", stage, begin, end, step, at, rep, grp, infoFields),
		m_rate(0), m_mode(mode), m_from(fromSubPop), m_to(toSubPop)
	{
		// when migrator is constructed from a pyMigrator, initial
		// rate is empty
		if (!rate.empty()) {
			DBG_FAILIF(!m_from.empty() && m_from.size() != rate.size(),
				ValueError, "Length of param fromSubPop must match rows of rate matrix.");

			DBG_FAILIF(!m_to.empty() && m_to.size() != rate[0].size(),
				ValueError, "Length of param toSubPop must match columns of rate matrix.");

			setRates(rate, mode);
		}
	};

	/// destructor
	virtual ~migrator()
	{
	};

	/// deep copy of a migrator
	virtual baseOperator * clone() const
	{
		return new migrator(*this);
	}


	/// return migration rate
	matrix rate()
	{
		return m_rate;
	}


	/// set migration rate
	/**
	   Format should be <tt>0-0 0-1 0-2, 1-0 1-1 1-2, 2-0, 2-1, 2-2</tt>.
	   For mode \c MigrByProbability or \c MigrByProportion, <tt>0-0,1-1,2-2</tt> will be set
	   automatically regardless of input.
	 */
	void setRates(const matrix & rate, int mode);

	/// apply the \c migrator
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c migrator
	virtual string __repr__()
	{
		return "<simuPOP::migrator>" ;
	}


protected:
	/// migration rate. its meaning is controled by m_mode
	matrix m_rate;

	/// asProbability (1), asProportion (2), or asCounts.
	int m_mode;

	/// from->to subPop index.
	/// default to 0 - rows of rate - 1, 0 - columns of rate - 1
	vectorvsp m_from;
	vectoru m_to;
};

/// a more flexible Python migrator
/**
   This migrator can be used in two ways
 \li define a function that accepts a generation number and returns a migration rate matrix.
   	This can be used in various migration rate cases.
 \li define a function that accepts individuals etc, and returns the new subpopulation ID.

   More specifically, \c func can be
 \li <tt>func(ind)</tt> when neither \c loci nor \c param is given.
 \li <tt>func(ind, genotype)</tt> when \c loci is given.
 \li <tt>func(ind, param)</tt> when \c param is given.
 \li <tt>func(ind, genotype, param)</tt> when both \c loci and \c param are given.

 */
class pyMigrator : public migrator
{

public:
	/// create a hybrid migrator
	/**
	 \param rateFunc a Python function that accepts a generation number,
	   	current subpopulation sizes, and returns a migration rate matrix.
	   	The migrator then migrate like a usual migrator.
	 \param indFunc a Python function that accepts an individual, optional
	   	genotypes and parameters, then returns a subpopulation ID. This
	   	method can be used to separate a population according to individual
	   	genotype.
	 \param stage default to \c PreMating
	 */
	pyMigrator(PyObject * rateFunc = NULL, PyObject * indFunc=NULL,
				int mode = MigrByProbability,
	           vectorvsp fromSubPop = vectorvsp(), vectoru toSubPop = vectoru(),
	           const vectoru & loci = vectoru(), PyObject * param = NULL,
	           int stage = PreMating, int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
	           int rep = REP_ALL, int grp = GRP_ALL, const vectorstr & infoFields = vectorstr())
		: migrator(matrix(), mode, fromSubPop, toSubPop, stage, begin, end, step, at, rep, grp, infoFields),
		m_rateFunc(rateFunc), m_indFunc(indFunc), m_loci(loci), m_param(param)
	{
		// carray of python list/typle
		DBG_FAILIF(rateFunc != NULL && !PyCallable_Check(rateFunc),
			ValueError, "Passed rateFunc is not a callable Python function.");
		DBG_FAILIF(indFunc != NULL && !PyCallable_Check(indFunc),
			ValueError, "Passed indFunc is not a callable Python function.");
		DBG_FAILIF(rateFunc == NULL && indFunc == NULL,
			ValueError, "Please specify either rateFunc or indFunc");
		DBG_FAILIF(rateFunc != NULL && indFunc != NULL,
			ValueError, "Please specify only one of rateFunc or indFunc");

		if (m_rateFunc != NULL)
			Py_INCREF(m_rateFunc);
		if (m_indFunc != NULL)
			Py_INCREF(m_indFunc);
		if (m_param != NULL)
			Py_INCREF(m_param);
	}


	/// destructor
	virtual ~pyMigrator()
	{
		if (m_rateFunc != NULL)
			Py_DECREF(m_rateFunc);
		if (m_indFunc != NULL)
			Py_DECREF(m_indFunc);
		if (m_param != NULL)
			Py_DECREF(m_param);
	}


	/// CPPONLY
	pyMigrator(const pyMigrator & rhs) : migrator(rhs),
		m_rateFunc(rhs.m_rateFunc),
		m_indFunc(rhs.m_indFunc),
		m_param(rhs.m_param)
	{
		if (m_rateFunc != NULL)
			Py_INCREF(m_rateFunc);
		if (m_indFunc != NULL)
			Py_INCREF(m_indFunc);
		if (m_param != NULL)
			Py_INCREF(m_param);
	}


	/// deep copy of a \c pyMigrator
	virtual baseOperator * clone() const
	{
		return new pyMigrator(*this);
	}


	/// apply a \c pyMigrator
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c pyMigrator
	virtual string __repr__()
	{
		return "<simuPOP::python migrator>" ;
	}


private:
	/// rateFunc
	PyObject * m_rateFunc;

	/// indFunc as in pyIndOperator
	PyObject * m_indFunc;

	/// loci to indFunc
	vectoru m_loci;

	/// parameters to indFunc
	PyObject * m_param;

};

/// split a subpopulation
/**
   <funcForm>SplitSubPop</funcForm>
 */
class splitSubPop : public baseOperator
{

public:
	/// split a subpopulation
	/**
	 \param which which subpopulation to split. If there is no subpopulation structure,
	   	use \c 0 as the first (and only) subpopulation.
	 \param sizes new subpopulation sizes. The sizes should be added up to the original
	   	subpopulation (subpopulation \c which) size.
	 \param proportions proportions of new subpopulations. Should be added up to \c 1.
	 \param subPopID new subpopulation IDs. Otherwise, the operator will automatically
	   	set new subpopulation IDs to new subpopulations.
	 \test src_splitSubPop.log Operator \c splitSubPop
	 */
	splitSubPop(UINT which = 0,  vectorlu sizes = vectorlu(), vectorf proportions = vectorf(),
	            vectoru subPopID = vectoru(),
	            bool randomize = true,
	            int stage = PreMating, int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
	            int rep = REP_ALL, int grp = GRP_ALL, const vectorstr & infoFields = vectorstr())
		: baseOperator("", "", stage, begin, end, step, at, rep, grp, infoFields),
		m_which(which), m_subPopSizes(sizes), m_proportions(proportions),
		m_subPopID(subPopID), m_randomize(randomize)
	{
		DBG_FAILIF(sizes.empty() && proportions.empty(), ValueError,
			"Please specify one of subPop and proportions.");
		DBG_FAILIF(!sizes.empty() && !proportions.empty(), ValueError,
			"Please specify only one of subPop and proportions.");
	}


	/// destructor
	virtual ~splitSubPop()
	{
	}


	/// deep copy of a \c splitSubPop operator
	virtual baseOperator * clone() const
	{
		return new splitSubPop(*this);
	}


	/// apply a \c splitSubPop operator
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c splitSubPop operator
	virtual string __repr__()
	{
		return "<simuPOP::split population>" ;
	}


private:
	/// which subpop to split
	UINT m_which;

	/// new subpopulation size
	vectorlu m_subPopSizes;

	/// new subpopulation proportions.
	vectorf m_proportions;

	/// subpopulation id, optional
	vectoru m_subPopID;

	/// random split
	/// randomize population before split.
	/// this is because some mating schemes generate
	/// individuals non-randomly, for example,
	/// put affected individuals at the beginning.
	bool m_randomize;
};

///  merge subpopulations
/**
   This operator merges subpopulations \c subPops to a
   single subpopulation. If \c subPops is ignored, all subpopulations will be merged.
   <funcForm>MergeSubPops</funcForm>
 */
class mergeSubPops : public baseOperator
{

public:
	/// merge subpopulations
	/**
	 \param subPops subpopulations to be merged. Default to all.
	 */
	mergeSubPops(vectoru subPops = vectoru(), bool removeEmptySubPops = false,
	             int stage = PreMating, int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
	             int rep = REP_ALL, int grp = GRP_ALL, const vectorstr & infoFields = vectorstr())
		: baseOperator("", "", stage, begin, end, step, at, rep, grp, infoFields),
		m_subPops(subPops), m_removeEmptySubPops(removeEmptySubPops)
	{
	}


	/// destructor
	virtual ~mergeSubPops()
	{
	}


	/// deep copy of a \c mergeSubPops operator
	virtual baseOperator * clone() const
	{
		return new mergeSubPops(*this);
	}


	/// apply a \c mergeSubPops operator
	virtual bool apply(population & pop)
	{
		pop.mergeSubPops(m_subPops, m_removeEmptySubPops);
		return true;
	}


	/// used by Python print function to print out the general information of the \c mergeSubPops operator
	virtual string __repr__()
	{
		return "<simuPOP::merge subpopulations>" ;
	}


private:
	///
	vectoru m_subPops;

	///
	bool m_removeEmptySubPops;
};


/// resize subpopulations
/**
   This operator resize subpopulations \c subPops to a
   another size. If \c subPops is ignored, all subpopulations will be resized.
   If the new size is smaller than the original one, the remaining individuals
   are discarded. If the new size if greater, individuals will be copied
   again if propagate is true, and be empty otherwise.
   <funcForm>ResizeSubPops</funcForm>
 */
class resizeSubPops : public baseOperator
{

public:
	/// resize subpopulations
	/**
	 \param newSizes of the specified (or all) subpopulations.
	 \param subPops subpopulations to be resized. Default to all.
	 \param propagate if true (default) and the new size if greater than
	   	the original size, individuals will be copied over.
	 */
	resizeSubPops(vectorlu newSizes = vectorlu(), vectoru subPops = vectoru(), bool propagate = true,
	              int stage = PreMating, int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
	              int rep = REP_ALL, int grp = GRP_ALL, const vectorstr & infoFields = vectorstr())
		: baseOperator("", "", stage, begin, end, step, at, rep, grp, infoFields),
		m_newSizes(newSizes), m_subPops(subPops), m_propagate(propagate)
	{
		DBG_FAILIF(!subPops.empty() && subPops.size() != newSizes.size(), ValueError,
			"Please specify new sizes for each specified subpopulation");
	}


	/// destructor
	virtual ~resizeSubPops()
	{
	}


	/// deep copy of a \c resizeSubPops operator
	virtual baseOperator * clone() const
	{
		return new resizeSubPops(*this);
	}


	/// apply a \c resizeSubPops operator
	virtual bool apply(population & pop);


	/// used by Python print function to print out the general information of the \c resizeSubPops operator
	virtual string __repr__()
	{
		return "<simuPOP::resize subpopulations>" ;
	}


private:
	///
	vectorlu m_newSizes;

	///
	vectoru m_subPops;

	///
	bool m_propagate;
};

}
#endif
