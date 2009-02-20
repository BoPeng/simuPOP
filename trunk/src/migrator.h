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


/** This operator migrates individuals from (virtual) subpopulations to other
 *  subpopulations, according to either pre-specified destination
 *  subpopulation stored in an information field, or randomly according to a
 *  migration matrix.
 *
 *  In the former case, values in a specified information field (default to
 *  \e migrate_to) are considered as destination subpopulation for each
 *  individual. If \e subPops is given, only individuals in specified (virtual)
 *  subpopulations will be migrated where others will stay in their original
 *  subpopulation. Negative values are not allowed in this information field
 *  because they do not represent a valid destination subpopulation ID.
 *
 *  In the latter case, a migration matrix is used to randomly assign
 *  destination subpoulations to each individual. The elements in this matrix
 *  can be probabilities to migrate, proportions of individuals to migrate, or
 *  exact number of individuals to migrate.
 *
 *  By default, the migration matrix should have \c m by \c m elements if there
 *  are \c m subpopulations. Element <tt>(i, j)</tt> in this matrix represents
 *  migration probability, rate or count from subpopulation \c i to \c j. If
 *  \e subPops (length \c m) and/or \e toSubPops (length \c n) are given,
 *  the matrix should have \c m by \c n elements, corresponding to specified
 *  source and destination subpopulations. Subpopulations in \e subPops can
 *  be virtual subpopulations, which makes it possible to migrate, for example,
 *  males and females at different rates from a subpopulation. If a
 *  subpopulation in \e toSubPops does not exist, it will be created. In case
 *  that all individuals from a subpopulation are migrated, the empty
 *  subpopulation will be kept.
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
 *  detination subpopulation of each individual to an information field. An
 *  \c RuntimeError will be raised if an individual is assigned to migrate
 *  more than once. This might happen if you are migrating from two overlapping
 *  virtual subpopulations.
 */
class migrator : public baseOperator
{
public:
	/** Create a migrator that moves individuals from source (virtual)
	 *  subpopulations \e subPops (default to migrate from all subpopulations)
	 *  to destination subpopulations \e toSubPops (default to all
	 *  subpopulations), according to existing values in an information field
	 *  \e infoFields[0], or randomly according to a migration matrix \e rate.
	 *  In the latter case, the size of the matrix should match the number of
	 *  source and destination subpopulations.
	 *
	 *  Depending on the value of parameter \e mode, elements in the migration
	 *  matrix (\e rate) are interpreted as either the probabilities to migrate
	 *  from source to destination subpopulations (\e mode = \c ByProbability),
	 *  proportions of individuals in the source (virtual) subpopulations to
	 *  the destination subpopulations (\e mode = \c ByProportion), numbers
	 *  of migrants in the source (virtual) subpopulations (\e mode
	 *  = \c ByCounts), or ignored completely (\e mode = \c ByIndInfo).
	 *  In the last case, parameter \e subPops is respected (only individuals
	 *  in specified (virtual) subpopulations will migrate) but \e toSubPops
	 *  is ignored.
	 *
	 *  This operator is by default applied pre-mating (parameter \e stage).
	 *  Please refer to operator \c baseOperator for a detailed explanation for
	 *  all parameters.
	 */
	migrator(const matrix & rate = matrix(), int mode = ByProbability, const uintList & toSubPops = uintList(),
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


	/** CPPONLY  set migration rate.
	 */
	void setRates(int mode, const subPopList & fromSubPops, const vectorlu & toSubPops);

	/// apply the migrator to populaiton \e pop.
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


/** Split a given list of subpopulations according to either sizes of the
 *  resulting subpopulations, proportion of individuals, or an information
 *  field. The resulting subpopulations will have the same name as the
 *  original subpopulation.
 *  <funcForm>SplitSubPops</funcForm>
 */
class splitSubPops : public baseOperator
{

public:
	/** Split a list of subpopulations \e subPops into finer subpopulations. A
	 *  single subpopulation is acceptable but virtual subpopulations are not
	 *  allowed. All subpopulations will be split if \e subPops is not specified.
	 *
	 *  The subpopulations can be split in three ways:
	 *  \li If parameter \e sizes is given, each subpopulation will be split
	 *  into subpopulations with given size. The \e sizes should add up to the
	 *  size of all orignal subpopulations.
	 *  \li If parameter \e proportions is given, each subpopulation will be
	 *  split into subpopulations with corresponding proportion of individuals.
	 *  \e proportions should add up to \c 1.
	 *  \li If an information field is given (parameter \e infoFields),
	 *  individuals having the same value at this information field will be
	 *  grouped into a subpopulation. The number of resulting subpopulations is
	 *  determined by the number of distinct values at this information field.
	 *
	 *  If parameter \c randomize is \c True (default), individuals will be
	 *  randomized before a subpopulation is split. This is designed to remove
	 *  artificial order of individuals introduced, for example, by some non-
	 *  random mating schemes. Note that, however, the original individual
	 *  order is not guaranteed even if this parameter is et to \c False.
	 *
	 *  Unless the last subpopulation is split, the indexes of existing
	 *  subpopulations will be changed. If a subpopulation has a name, this
	 *  name will become the name for all subpopulations separated from this
	 *  subpopulation.
	 *
	 *  This operator is by default applied pre-mating (parameter \e stage).
	 *  Please refer to operator \c baseOperator for a detailed explanation for
	 *  all parameters.
	 *
	 *  \note Unlike operator \c migrator, this operator does not require an
	 *  information field such as \c migrate_to.
	 */
	splitSubPops(const subPopList & subPops = subPopList(), const vectorlu & sizes = vectorlu(),
		const vectorf & proportions = vectorf(), bool randomize = true,
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const vectorstr & infoFields = vectorstr())
		: baseOperator("", stage, begin, end, step, at, rep, subPops, infoFields),
		m_subPopSizes(sizes), m_proportions(proportions), m_randomize(randomize)
	{
		for (size_t i = 0; i < subPops.size(); ++i) {
			DBG_FAILIF(subPops[i].isVirtual(), ValueError,
				"Virtual subpopulations are not supported in operator splitSubPops");
		}

		DBG_FAILIF(sizes.empty() + proportions.empty() + infoFields.empty() != 2, ValueError,
			"Please specify one and only one of sizes, proportions and infoFields.");
	}


	/// destructor
	virtual ~splitSubPops()
	{
	}


	/// deep copy of a \c splitSubPops operator
	virtual baseOperator * clone() const
	{
		return new splitSubPops(*this);
	}


	/// apply a \c splitSubPops operator
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c splitSubPops operator
	virtual string __repr__()
	{
		return "<simuPOP::split subpopulation>" ;
	}


private:
	/// new subpopulation size
	vectorlu m_subPopSizes;

	/// new subpopulation proportions.
	vectorf m_proportions;

	/// random split
	bool m_randomize;
};


/** This operator merges subpopulations \e subPops to a single subpopulation.
 *  If \c subPops is ignored, all subpopulations will be merged. Virtual
 *  subpopulations are not allowed in \e subPops.
 *  <funcForm>MergeSubPops</funcForm>
 */
class mergeSubPops : public baseOperator
{

public:
	/** Create an operator that merges subpopulations \e subPops to a single
	 *  subpopulation. If \e subPops is not given, all subpopulations will be
	 *  merged. The merged subpopulation will take the name of the first
	 *  subpopulation being merged.
	 *
	 *  This operator is by default applied pre-mating (parameter \e stage).
	 *  Please refer to operator \c baseOperator for a detailed explanation for
	 *  all parameters.
	 */
	mergeSubPops(const subPopList & subPops = subPopList(),
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const vectorstr & infoFields = vectorstr())
		: baseOperator("", stage, begin, end, step, at, rep, subPops, infoFields)
	{
		for (size_t i = 0; i < subPops.size(); ++i) {
			DBG_FAILIF(subPops[i].isVirtual(), ValueError,
				"Virtual subpopulations are not supported in operator mergeSubPops");
		}
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
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c mergeSubPops operator
	virtual string __repr__()
	{
		return "<simuPOP::merge subpopulations>" ;
	}


};


/** This operator resizes subpopulations to specified sizes. Individuals are
 *  added or removed depending on the new subpopulation sizes.
 *  <funcForm>ResizeSubPops</funcForm>
 */
class resizeSubPops : public baseOperator
{

public:
	/** Resize given subpopulations \e subPops to new sizes \e size, or sizes
	 *  proportional to original sizes (parameter \e proportions). All
	 *  subpopulations will be resized if \e subPops is not specified. If the
	 *  new size of a subpopulation is smaller than its original size, extra
	 *  individuals will be removed. If the new size is larger, new individuals
	 *  with empty genotype will be inserted, unless parameter \e propagate is
	 *  set to \c True (default). In this case, existing individuals will be
	 *  copied sequentially, and repeatedly if needed.
	 *
	 *  This operator is by default applied pre-mating (parameter \e stage).
	 *  Please refer to operator \c baseOperator for a detailed explanation for
	 *  all parameters.
	 */
	resizeSubPops(const subPopList & subPops = subPopList(),
		const vectorlu & sizes = vectorlu(), const vectorf & proportions = vectorf(),
		bool propagate = true, int stage = PreMating,
		int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const vectorstr & infoFields = vectorstr())
		: baseOperator("", stage, begin, end, step, at, rep, subPops, infoFields),
		m_sizes(sizes), m_proportions(proportions), m_propagate(propagate)
	{
		for (size_t i = 0; i < subPops.size(); ++i) {
			DBG_FAILIF(subPops[i].isVirtual(), ValueError,
				"Virtual subpopulations are not supported in operator splitSubPops");
		}

		DBG_FAILIF(m_sizes.empty() + m_proportions.empty() != 1, ValueError,
			"Please specify one and only one of parameter sizes and proportions");

		DBG_FAILIF(!subPops.empty() &&
			((subPops.size() != m_sizes.size()) && (subPops.size() != m_proportions.size())),
			ValueError,
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
	vectorlu m_sizes;

	vectorf m_proportions;
	///
	bool m_propagate;
};

}
#endif
