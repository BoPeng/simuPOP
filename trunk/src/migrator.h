/**
 *  $File: migrator.h $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MIGRATOR_H
#define _MIGRATOR_H
/**
   \file
   \brief head file of class Migrator:public BaseOperator
 */
#include "operator.h"
#include <list>
using std::list;

#include "boost_pch.hpp"


#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
using std::sort;

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
 *  detination subpopulation of each individual to an information field.
 *  Unexpected results may happen if individuals migrate from overlapping
 *  virtual subpopulations.
 */
class Migrator : public BaseOperator
{
public:
	/** Create a Migrator that moves individuals from source (virtual)
	 *  subpopulations \e subPops (default to migrate from all subpopulations)
	 *  to destination subpopulations \e toSubPops (default to all
	 *  subpopulations), according to existing values in an information field
	 *  \e infoFields[0], or randomly according to a migration matrix \e rate.
	 *  In the latter case, the size of the matrix should match the number of
	 *  source and destination subpopulations.
	 *
	 *  Depending on the value of parameter \e mode, elements in the migration
	 *  matrix (\e rate) are interpreted as either the probabilities to migrate
	 *  from source to destination subpopulations (\e mode = \c BY_PROBABILITY),
	 *  proportions of individuals in the source (virtual) subpopulations to
	 *  the destination subpopulations (\e mode = \c BY_PROPORTION), numbers
	 *  of migrants in the source (virtual) subpopulations (\e mode
	 *  = \c BY_COUNTS), or ignored completely (\e mode = \c BY_IND_INFO).
	 *  In the last case, parameter \e subPops is respected (only individuals
	 *  in specified (virtual) subpopulations will migrate) but \e toSubPops
	 *  is ignored.
	 *
	 *  Please refer to operator \c BaseOperator for a detailed explanation for
	 *  all parameters.
	 */
	Migrator(const floatMatrix & rate = floatMatrix(), int mode = BY_PROBABILITY,
		const uintList & toSubPops = uintList(),
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "migrate_to"));

	/// destructor
	virtual ~Migrator()
	{
	};

	/// HIDDEN Deep copy of a Migrator
	virtual BaseOperator * clone() const
	{
		return new Migrator(*this);
	}


	/// HIDDEN apply the Migrator to populaiton \e pop.
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const;

protected:
	/// migration rate. its meaning is controled by m_mode
	const matrixf m_rate;

	/// asProbability (1), asProportion (2), or asCounts.
	const int m_mode;

	/// from->to subPop index.
	/// default to 0 - rows of rate - 1, 0 - columns of rate - 1
	const uintList m_to;
};


/** This operator migrates individuals between all available or specified 
 *  subpopulations, according to a backward migration matrix. It differs from
 *  \c Migrator in how migration matrixes are interpreted. Due to the limit
 *  of this model, this operator does not support migration by information
 *  field, migration by count (\e mode = \c BY_COUNT), migration from virtual
 *  subpopulations, migration between different number of subpopulations,
 *  and the creation of new subpopulation, as operator \c Migrator provides.
 *
 *  In contrast to a forward migration matrix where $m_{ij}$ is considered
 *  the probability (proportion or count) of individuals migrating from subpopulation
 *  \c i to \c j, elements in a reverse migration matrix $m_{ij}$ is considered
 *  the probability (proportion or count) of individuals migrating from subpopulation
 *  \c j to \c i, namely the probability (proportion or count) of individuals
 *  originats from subpopulation \c j. 
 *
 *  If migration is applied by probability, the row of the migration matrix
 *  corresponding to a destination subpopulation is intepreted as probabilities to
 *  orignate from each source subpopulation. Each individual's source
 *  subpopulation is assigned randomly according to these probabilities. Note
 *  that the probability of originating from the present subpopulation is
 *  automatically calculated so the corresponding matrix elements are ignored.
 *
 *  If migration is applied by proportion, the row of the migration matrix
 *  corresponding to a destination subpopulation is intepreted as proportions
 *  to originate from each source subpopulation. The number of migrants from each
 *  source subpopulation is determined before random indidividuals are
 *  chosen to migrate.
 *
 *  Unlike the forward migration matrix that describes how migration should
 *  be performed, the backward migration matrix describes the result of
 *  migration. The underlying forward migration matrix is calculated at
 *  each generation and is in theory not the same across generations.
 * 
 *  This operator calculates the corresponding forward migration matrix
 *  from backward matrix and current population size. This process is not
 *  always feasible so an error will raise if no valid ending population
 *  size or forward migration matrix could be determined. Please refer to 
 *  the simuPOP user's guide for an explanation of the theory behind forward
 *  and backward migration matrices.
 */
class BackwardMigrator : public BaseOperator
{
public:
	/** Create a BackwardMigrator that moves individuals between \e subPop
	 *  subpopulations randomly according to a backward migration matrix \e rate.
	 *  The size of the matrix should match the number of subpopulations.
	 *
	 *  Depending on the value of parameter \e mode, elements in the migration
	 *  matrix (\e rate) are interpreted as either the probabilities to originate
	 *  from source subpopulations (\e mode = \c BY_PROBABILITY) or proportions of
	 *  individuals originate from the source (virtual) subpopulations (\e mode
	 *  = \c BY_PROPORTION). Migration by count is not supported by this operator.
	 *
	 *  Please refer to operator \c BaseOperator for a detailed explanation for
	 *  all parameters.
	 */
	BackwardMigrator(const floatMatrix & rate = floatMatrix(), int mode = BY_PROBABILITY,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "migrate_to"));

	/// destructor
	virtual ~BackwardMigrator()
	{
	};

	/// HIDDEN Deep copy of a Migrator
	virtual BaseOperator * clone() const
	{
		return new BackwardMigrator(*this);
	}


	/// HIDDEN apply the Migrator to populaiton \e pop.
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const;

protected:
	/// migration rate. its meaning is controled by m_mode
	const matrixf m_rate;

	boost::numeric::ublas::matrix<double> m_inverse_rate;

	bool m_symmetric_matrix;

	/// asProbability (1), asProportion (2),
	const int m_mode;
};




/** Split a given list of subpopulations according to either sizes of the
 *  resulting subpopulations, proportion of individuals, or an information
 *  field. The resulting subpopulations will have the same name as the
 *  original subpopulation.
 */
class SplitSubPops : public BaseOperator
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
	 *  artificial order of individuals introduced by, for example, some non-
	 *  random mating schemes. Note that, however, the original individual
	 *  order is not guaranteed even if this parameter is set to \c False.
	 *
	 *  Unless the last subpopulation is split, the indexes of existing
	 *  subpopulations will be changed. If a subpopulation has a name, this
	 *  name will become the name for all subpopulations separated from this
	 *  subpopulation. Optionally, you can assign names to the new
	 *  subpopulations using a list of names specified in parameter \e names.
	 *  Because the same set of names will be used for all subpopulations,
	 *  this parameter is not recommended when multiple subpopulations are
	 *  split.
	 *
	 *  Please refer to operator \c BaseOperator for a detailed explanation for
	 *  all parameters.
	 *
	 *  \note Unlike operator \c Migrator, this operator does not require an
	 *  information field such as \c migrate_to.
	 */
	SplitSubPops(const subPopList & subPops = subPopList(), const vectoru & sizes = vectoru(),
		const vectorf & proportions = vectorf(), const stringList names = vectorstr(), bool randomize = true,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const stringList & infoFields = vectorstr())
		: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_subPopSizes(sizes), m_proportions(proportions), m_names(names.elems()), m_randomize(randomize)
	{
		for (size_t i = 0; i < subPops.size(); ++i) {
			DBG_FAILIF(subPops[i].isVirtual(), ValueError,
				"Virtual subpopulations are not supported in operator SplitSubPops");
		}

		DBG_ASSERT(m_proportions.empty() || fcmp_eq(accumulate(m_proportions.begin(), m_proportions.end(), 0.), 1.),
			ValueError, "Proportions should add up to one.");
		DBG_FAILIF(sizes.empty() + proportions.empty() + infoFields.elems().empty() != 2, ValueError,
			"Please specify one and only one of sizes, proportions and infoFields.");
	}


	/// destructor
	virtual ~SplitSubPops()
	{
	}


	/// HIDDEN Deep copy of a \c SplitSubPops operator
	virtual BaseOperator * clone() const
	{
		return new SplitSubPops(*this);
	}


	/// HIDDEN apply a \c SplitSubPops operator
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.split subpopulation>" ;
	}


private:
	/// new subpopulation size
	vectoru m_subPopSizes;

	/// new subpopulation proportions.
	vectorf m_proportions;

	vectorstr m_names;

	/// random split
	bool m_randomize;
};


/** This operator merges subpopulations \e subPops to a single subpopulation.
 *  If \c subPops is ignored, all subpopulations will be merged. Virtual
 *  subpopulations are not allowed in \e subPops.
 */
class MergeSubPops : public BaseOperator
{

public:
	/** Create an operator that merges subpopulations \e subPops to a single
	 *  subpopulation. If \e subPops is not given, all subpopulations will be
	 *  merged. The merged subpopulation will take the name of the first
	 *  subpopulation being merged unless a new \e name is given.
	 *
	 *  Please refer to operator \c BaseOperator for a detailed explanation for
	 *  all parameters.
	 */
	MergeSubPops(const subPopList & subPops = subPopList(), const string & name = string(),
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const stringList & infoFields = vectorstr())
		: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_name(name)
	{
		for (size_t i = 0; i < subPops.size(); ++i) {
			DBG_FAILIF(subPops[i].isVirtual(), ValueError,
				"Virtual subpopulations are not supported in operator MergeSubPops");
		}
	}


	/// destructor
	virtual ~MergeSubPops()
	{
	}


	/// HIDDEN Deep copy of a \c MergeSubPops operator
	virtual BaseOperator * clone() const
	{
		return new MergeSubPops(*this);
	}


	/// HIDDEN apply a \c MergeSubPops operator
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.merge subpopulations>" ;
	}


private:
	string m_name;
};


/** This operator resizes subpopulations to specified sizes. individuals are
 *  added or removed depending on the new subpopulation sizes.
 */
class ResizeSubPops : public BaseOperator
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
	 *  Please refer to operator \c BaseOperator for a detailed explanation for
	 *  all parameters.
	 */
	ResizeSubPops(const subPopList & subPops = subPopList(),
		const vectoru & sizes = vectoru(), const vectorf & proportions = vectorf(),
		bool propagate = true,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const stringList & infoFields = vectorstr())
		: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_sizes(sizes), m_proportions(proportions), m_propagate(propagate)
	{
		for (size_t i = 0; i < subPops.size(); ++i) {
			DBG_FAILIF(subPops[i].isVirtual(), ValueError,
				"Virtual subpopulations are not supported in operator SplitSubPops");
		}

		DBG_FAILIF(m_sizes.empty() + m_proportions.empty() != 1, ValueError,
			"Please specify one and only one of parameter sizes and proportions");

		DBG_FAILIF(!subPops.empty() &&
			((subPops.size() != m_sizes.size()) && (subPops.size() != m_proportions.size())),
			ValueError,
			"Please specify new sizes for each specified subpopulation");
	}


	/// destructor
	virtual ~ResizeSubPops()
	{
	}


	/// HIDDEN Deep copy of a \c ResizeSubPops operator
	virtual BaseOperator * clone() const
	{
		return new ResizeSubPops(*this);
	}


	/// HIDDEN apply a \c ResizeSubPops operator
	virtual bool apply(Population & pop) const;


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.resize subpopulations>" ;
	}


private:
	///
	vectoru m_sizes;

	vectorf m_proportions;
	///
	bool m_propagate;
};

}
#endif
