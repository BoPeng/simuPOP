/**
 *  $File: qtrait.h $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
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

#ifndef _QTRAIT_H
#define _QTRAIT_H
/**
   \file
   \brief head file of class selector:public baseOperator
 */
#include "utility.h"
#include "operator.h"

#include "boost/tuple/tuple.hpp"
#include <numeric>
using std::min;

namespace simuPOP {

// ///////////////////////// Quantitative trait ////////////////////////

/// base class of quantitative trait
/**
   Quantitative trait is the measure of certain phenotype for given genotype.
   Quantitative trait is similar to penetrance in that the consequence of
   penetrance is binary: affected or unaffected; while it is continuous for
   quantitative trait. \n

   In simuPOP, different operators or functions were implemented to calculate
   quantitative traits for each individual and store the values in the information
   fields specified by the user (default to \c qtrait). The quantitative trait operators
   also accept the \c ancestralGen parameter to control the number of generations
   for which the \c qtrait information field will be set.
 */
class quanTrait : public baseOperator
{
public:
	/// create a quantitative trait operator
	quanTrait(int ancGen = -1,  int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("qtrait"))
		: baseOperator("", stage, begin, end, step, at, reps, subPops, infoFields),
		m_ancGen(ancGen)
	{
	}


	/// destructor
	virtual ~quanTrait()
	{
	}


	/// deep copy of a quantitative trait operator
	virtual baseOperator * clone() const
	{
		return new quanTrait(*this);
	}


	/** CPPONLY
	 *  calculate/return quantitative trait etc.
	 */
	virtual double qtrait(individual *)
	{
		///
		throw ValueError("This quantitative trait calculator is not supposed to be called directly");
		return 1.;
	}


	/// set \c qtrait to all individual
	bool apply(population & pop);

	/// HIDDEN
	string description()
	{
		return "<simuPOP.qtrait::quantitative trait>" ;
	}


private:
	/// how to handle ancestral gen
	int m_ancGen;

};

/// quantitative trait according to genotype at one locus
/**
   Assign quantitative trait using a table with keys 'X-Y' where X and Y are allele
   numbers. If parameter \c sigma is not zero, the return value is the sum of the
   trait plus \f$ N\left(0,\sigma^{2}\right) \f$. This random part is usually considered
   as the environmental factor of the trait.
   <funcForm>MapQuanTrait</funcForm>
 */
class mapQuanTrait : public quanTrait
{
public:
	/// create a map quantitative trait operator
	/**
	   \param locus the locus index. The quantitative trait is determined by genotype at this locus.
	   \param loci an array of locus indexes. The quantitative trait is determined by genotypes at these loci.
	   \param qtrait a dictionary of quantitative traits. The genotype must be in the
	    form of 'a-b'. This is the mean	of the quantitative trait. The actual trait
	    value will be \f$ N\left(mean,\sigma^{2}\right) \f$. For multiple loci, the form is
	    'a-b|c-d|e-f' etc.
	   \param sigma standard deviation of the environmental factor \f$ N\left(0,\sigma^{2}\right) \f$.
	   \param phase if \c True, a/b and b/a will have different quantitative trait values.
	    Default to \c False.
	   \param output and other parameters please refer to help (<tt>baseOperator.__init__</tt>)
	 */
	mapQuanTrait(const uintList & loci, const strDict & qtrait, double sigma = 0, bool phase = false,
		int ancGen = -1,
		int stage = PostMating, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("qtrait")) :
		quanTrait(ancGen, stage, begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci.elems()), m_dict(qtrait), m_sigma(sigma), m_phase(phase)
	{
	};

	virtual ~mapQuanTrait()
	{
	}


	/// deep copy of a map quantitative trait operator
	virtual baseOperator * clone() const
	{
		return new mapQuanTrait(*this);
	}


	/** CPPONLY
	 *  currently assuming diploid
	 */
	virtual double qtrait(individual * ind);

	/// HIDDEN
	string description()
	{
		return "<simuPOP.qtrait::map quantitative trait>" ;
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

/// multiple allele quantitative trait (quantitative trait according to disease or wildtype alleles)
/**
   This is called 'multiple-allele' quantitative trait. It separates alleles into
   two groups: wildtype and diseased alleles. Wildtype alleles are specified by parameter
   \c wildtype and any other alleles are considered as diseased alleles.
   \c maQuanTrait accepts an array of fitness. Quantitative trait is then set for any given
   genotype. A standard normal distribution \f$ N\left(0,\sigma^{2}\right) \f$ will
   be added to the returned trait value.
   <funcForm>MaQuanTrait</funcForm>
 */
class maQuanTrait : public quanTrait
{
public:
	/// create a multiple allele quantitative trait operator
	/**
	   \param qtrait an array of quantitative traits of AA, Aa, aa. A is the wildtype group
	   \param sigma an array of standard deviations for each of the trait genotype (AA, Aa, aa)
	   \param wildtype an array of alleles in the wildtype group. Any other alleles will be
	    considered as diseased alleles. Default to <tt>[0]</tt>.
	   \param output and other parameters please refer to help(<tt>baseOperator.__init__</tt>)

	   Please refer to \c quanTrait for other parameter descriptions.
	 */
	maQuanTrait(const uintList & loci, const vectorf & qtrait, const uintList & wildtype,
		const floatList & sigma = vectorf(), int ancGen = -1,
		int stage = PostMating, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("qtrait"));

	/// destructor
	virtual ~maQuanTrait()
	{
	}


	/// deep copy of a multiple allele quantitative trait
	virtual baseOperator * clone() const
	{
		return new maQuanTrait(*this);
	}


	/** CPPONLY
	 * currently assuming diploid
	 */
	virtual double qtrait(individual * ind);

	/// HIDDEN
	string description()
	{
		return "<simuPOP.qtrait::multiple-alleles qtrait>" ;
	}


private:
	/// one locus
	vectoru m_loci;

	/// qtrait for each genotype
	vectorf m_qtrait;

	///
	vectorf m_sigma;

	///
	vectoru m_wildtype;
};

/// quantitative trait according to genotypes from a multiple loci multiplicative model
/**
   Operator \c mlQuanTrait is a 'multiple-locus' quantitative trait calculator. It accepts a list
   of quantitative traits and combine them according to the \c mode parameter, which takes
   one of the following values
   \li \c Multiplicative: the mean of the quantitative trait is calculated as
   \f$ f=\prod f_{i} \f$.
   \li \c Additive: the mean of the quantitative trait is calculated as
   \f$ f=\sum f_{i} \f$.

   Note that all \f$ \sigma_{i} \f$ (for \f$ f_{i} \f$) and \f$ \sigma \f$ (for \f$ f \f$)
   will be considered. I.e, the trait value should be
   \f[ f=\sum_{i}\left(f_{i}+N\left(0,\sigma_{i}^{2}\right)\right)+\sigma^{2} \f]
   for \c Additive case. If this is not desired, you can set some of the \f$ \sigma \f$ to zero.

   <funcForm>MlQuanTrait</funcForm>
 */
class mlQuanTrait : public quanTrait
{

public:
	/// create a multiple locus quantitative trait operator
	/**
	   \param qtraits a list of quantitative traits
	   \param mode can be one of \c Multiplicative and \c Additive

	   Please refer to \c quanTrait for other parameter descriptions.
	 */
	mlQuanTrait(const opList & qtraits, int mode = Multiplicative,
		double sigma = 0, int ancGen = -1,
		int stage = PostMating, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("qtrait")) :
		quanTrait(ancGen, stage, begin, end, step, at, reps, subPops, infoFields),
		m_qtraits(qtraits), m_sigma(sigma), m_mode(mode)
	{
		DBG_FAILIF(qtraits.empty(), ValueError, "Please specify at least one selector.");
	};

	virtual ~mlQuanTrait()
	{
	}


	/// deep copy of a multiple loci quantitative trait operator
	virtual baseOperator * clone() const
	{
		throw ValueError("Multi-loci selector can not be nested.");
	}


	/// CPPONLY
	/// currently assuming diploid
	virtual double qtrait(individual * ind);

	/// HIDDEN
	string description()
	{
		return "<simuPOP.qtrait::multiple-loci qtrait>" ;
	}


private:
	/// a list of qtraits
	opList m_qtraits;

	///
	double m_sigma;

	/// mode
	int m_mode;
};

/// quantitative trait using a user provided function
/**
   For each individual, a user provided function is used to calculate quantitative trait.
   <funcForm>PyQuanTrait</funcForm>
 */
class pyQuanTrait : public quanTrait
{
public:
	/// create a Python quantitative trait operator
	/**
	   \param loci The genotypes at these loci will be
	    passed to \c func.
	   \param func a Python function that accepts genotypes at specified loci
	    and returns the quantitative trait value.
	   \param output and other parameters please refer to help(<tt>baseOperator.__init__</tt>)

	   Please refer to \c quanTrait for other parameter descriptions.
	 */
	// provide locus and qtrait for 11, 12, 13 (in the form of dictionary)
	pyQuanTrait(const uintList & loci, PyObject * func, int ancGen = -1,
		int stage = PostMating, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("qtrait")) :
		quanTrait(ancGen, stage, begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci.elems()), m_func(func), m_alleles(0), m_len(0), m_numArray(NULL)
	{
		if (!m_func.isValid())
			throw ValueError("Passed variable is not a callable python function.");

		DBG_FAILIF(m_loci.empty(), ValueError,
			"Please specify susceptibility loci");
	};


	/// CPPONLY
	pyQuanTrait(const pyQuanTrait & rhs) :
		quanTrait(rhs),
		m_loci(rhs.m_loci),
		m_func(rhs.m_func),
		m_alleles(rhs.m_alleles),
		m_len(rhs.m_len),
		m_numArray(NULL)
	{
	}


	/// deep copy of a Python quantitative trait operator
	virtual baseOperator * clone() const
	{
		return new pyQuanTrait(*this);
	}


	/** CPPONLY
	 *  currently assuming diploid
	 */
	virtual double qtrait(individual * ind);

	/// HIDDEN
	string description()
	{
		return "<simuPOP.qtrait::python qtrait>" ;
	}


private:
	/// susceptibility loci
	vectoru m_loci;

	/// user supplied python function
	pyFunc m_func;

	/// copy of alleles of each individual a time.
	vectora m_alleles;

	/// length of m_alleles
	int m_len;

	/// the object that passed to func
	PyObject * m_numArray;

};

}
#endif
