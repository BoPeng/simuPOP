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


/** This operator migrate individuals from (virtual) subpopulations to other
 *  subpopulations according to a migration matrix. The elements in this matrix
 *  can be probabilities to migrate, proportions of individuals to migrate,
 *  or exact number of individuals to migrate. An information field (default
 *  to \c migrate_to) is required for this operator.
 *
 *  By default, the migration matrix should have \c m by \c m elements if there
 *  are \c m subpopulations. Element <tt>(i, j)</tt> in this matrix represents
 *  migration probability, rate or count from subpopulation \c i to \c j. If
 *  \e subPops (length \c m) and/or \e toSubPops (length \c n) are given,
 *  the matrix should have \c m by \c n elements, corresponding to specified 
 *  source and destination subpopulations. Subpopulations in \e fromSubPops can
 *  be virtual subpopulations, which makes it possible to migrate, for example,
 *  males and females at different rates from a subpopulation.
 *  
 *  If migration is applied by probability, the row of the migration matrix
 *  corresponding to a source subpopulation is intepreted as probabilities to
 *  migrate to each destination subpopulation. Each individual's detination
 *  subpopulation is assigned randomly according to these probabilities. Note
 *  that the probability of staying at the present subpopulation is
 *  automatically calculated so the corresponding matrix elements are ignored.
 *
 *  If migration is applied by proportion, the row of the migration matrix
 *  corresponding to a source subpopulation is intepreted as proportions to
 *  migrate to each destination subpopulation. The number of migrants to each
 *  destination subpopulation is determined before random indidividuals are
 *  chosen to migrate.
 *
 *  If migration is applied by counts, the row of the migration matrix
 *  corresponding to a source subpopulation is intepreted as number of
 *  individuals to migrate to each detination subpopulation. The migrants are
 *  chosen randomly.
 *  
 *  This operator goes through all source (virtual) subpopulations and assign
 *  detination subpopulation of each individual to an information field. If
 *  overlapped or duplicated (virtual) source subpopulations are given, latter
 *  destination assignments will override previous assignments. After all
 *  source subpopulations are handled, individuals will migrate to their
 *  assigned destination subpopulation altogether.
 */
class migrator : public baseOperator
{
public:
	/** Create a migrator that moves individuals from source (virtual)
	 *  subpopulations \e subPops (default to migrate from all subpopulations)
	 *  to destination subpopulations \e toSubPops (default to all
	 *  subpopulations), according to a migration matrix \e rate. The size of
	 *  this matrix should match the number of source and destination (virtual)
	 *  subpopulations. This migrator goes through all source subpoulations,
	 *  assign each individual's destination subpopulation to an information
	 *  field \e infoFields[0] (default to \c migrate_to), and migrate
	 *  individuals to their destination subpopulations altogether.
	 *
	 *  Depending on the value of parameter \e mode, elements in the migration
	 *  matrix (\e rate) are interpreted as either the probabilities to migrate
	 *  from source to destination subpopulations (\e mode = \c ByProbability),
	 *  proportions of individuals in the source (virtual) subpopulations to
	 *  the destination subpopulations (\e mode = \c ByProportion), or numbers
	 *  of migrants in the source (virtual) subpopulations (\e mode
	 *  = \c ByCounts).
	 */
	migrator(const matrix & rate, int mode = ByProbability, const uintList & toSubPops = uintList(),
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(),
		const vectorstr & infoFields = vectorstr(1, "migrate_to"));

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
	   For mode \c ByProbability or \c ByProportion, <tt>0-0,1-1,2-2</tt> will be set
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
	vectorlu m_to;
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
	pyMigrator(PyObject * rateFunc = NULL, PyObject * indFunc = NULL,
		int mode = ByProbability, const uintList & toSubPops = uintList(),
		const vectoru & loci = vectoru(), PyObject * param = NULL,
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr(1, "migrate_to"))
		: migrator(matrix(), mode, toSubPops, stage, begin, end, step, at, rep, subPops, infoFields),
		m_rateFunc(rateFunc), m_indFunc(indFunc), m_loci(loci), m_param(param)
	{
		DBG_FAILIF(!m_rateFunc.isValid() && !m_indFunc.isValid(),
			ValueError, "Please specify either rateFunc or indFunc");
		DBG_FAILIF(m_rateFunc.isValid() && m_indFunc.isValid(),
			ValueError, "Please specify only one of rateFunc or indFunc");
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
	pyFunc m_rateFunc;

	/// indFunc as in pyIndOperator
	pyFunc m_indFunc;

	/// loci to indFunc
	vectoru m_loci;

	/// parameters to indFunc
	pyObject m_param;

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
	   Split a subpopulation by sizes or proportions. Individuals are randomly (by default)
	   assigned to the resulting subpopulations. Because mating schemes may introduce
	   certain order to individuals, randomization ensures that split subpopulations have
	   roughly even distribution of genotypes.

	   \param which which subpopulation to split. If there is no subpopulation structure,
	    use \c 0 as the first (and only) subpopulation.
	   \param sizes new subpopulation sizes. The sizes should be added up to the original
	    subpopulation (subpopulation \c which) size.
	   \param proportions proportions of new subpopulations. Should be added up to \c 1.
	   \param randomize Whether or not randomize individuals before population split. Default
	    to True.
	 */
	splitSubPop(UINT which = 0,  vectorlu sizes = vectorlu(), vectorf proportions = vectorf(),
		bool randomize = true,
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr(1, "migrate_to"))
		: baseOperator("", stage, begin, end, step, at, rep, subPops, infoFields),
		m_which(which), m_subPopSizes(sizes), m_proportions(proportions),
		m_randomize(randomize)
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
	mergeSubPops(const subPopList & subPops = subPopList(),
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const vectorstr & infoFields = vectorstr())
		: baseOperator("", stage, begin, end, step, at, rep, subPops, infoFields)
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
		subPopList sp = applicableSubPops();
		vectoru subPops(sp.size());

		for (size_t i = 0; i < sp.size(); ++i)
			subPops[i] = sp[i].subPop();
		pop.mergeSubPops(subPops);
		return true;
	}


	/// used by Python print function to print out the general information of the \c mergeSubPops operator
	virtual string __repr__()
	{
		return "<simuPOP::merge subpopulations>" ;
	}


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
	   \param subPop subpopulations to be resized. Default to all.
	   \param propagate if true (default) and the new size if greater than
	    the original size, individuals will be copied over.
	 */
	resizeSubPops(vectorlu newSizes = vectorlu(), bool propagate = true,
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: baseOperator("", stage, begin, end, step, at, rep, subPops, infoFields),
		m_newSizes(newSizes), m_propagate(propagate)
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
	bool m_propagate;
};

}
#endif
