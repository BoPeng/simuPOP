/**
 *  $File: virtualSubPop.h $
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

#ifndef _VIRTUALSUBPOP_H
#define _VIRTUALSUBPOP_H
/**
   \file
   \brief head file of class mating and its subclasses
 */
#include "utility.h"
#include "simuPOP_cfg.h"

#include <cmath>

namespace simuPOP {

class Individual;
class Population;

/**
 * A class to specify virtual subpopulation, which is composed of a subpopulation
 * ID and a virtual subpopulation ID.
 *
 */
class vspID
{
public:
	/** Create a subpopulation id. Accept id as well as names. */
	vspID(PyObject * id);

	/// CPPONLY
	vspID(const vectori & subPop, bool allAvailSP = false, bool allAvailVSP = false,
		const string & spName = string(), const string & vspName = string()) :
		m_spName(spName), m_vspName(vspName), m_allAvailSP(allAvailSP), m_allAvailVSP(allAvailVSP)
	{
		DBG_FAILIF(subPop.size() > 2, ValueError,
			"VSP should be specified as a subPop and virtualSubPop ID pair");
		m_subPop = subPop.size() > 0 && subPop[0] >= 0 ? subPop[0] : InvalidValue;
		m_virtualSubPop = subPop.size() > 1 && subPop[1] >= 0 ? subPop[1] : InvalidValue;
	}


	/// CPPONLY
	vspID(size_t subPop = InvalidValue, size_t virtualSubPop = InvalidValue,
		bool allAvailSP = false, bool allAvailVSP = false,
		const string & spName = string(), const string & vspName = string())
		: m_subPop(subPop), m_virtualSubPop(virtualSubPop),
		m_spName(spName), m_vspName(vspName), m_allAvailSP(allAvailSP), m_allAvailVSP(allAvailVSP)
	{
	}


	~vspID()
	{
	}


	/// CPPONLY
	bool operator==(const vspID & rhs) const
	{
		return m_subPop == rhs.m_subPop && m_virtualSubPop == rhs.m_virtualSubPop;
	}


	/// CPPONLY
	size_t subPop() const
	{
		DBG_FAILIF(!m_spName.empty(), ValueError, "Unresolved population name.");
		return m_subPop;
	}


	/// CPPONLY
	size_t virtualSubPop() const
	{
		DBG_FAILIF(!m_vspName.empty(), ValueError, "Unresolved virtual subpopulation name.");
		return m_virtualSubPop;
	}


	/// CPPONLY
	bool valid() const
	{
		return m_subPop != InvalidValue;
	}


	/// CPPONLY
	bool isVirtual() const
	{
		return m_virtualSubPop != InvalidValue;
	}


	/// CPPONLY
	bool allAvailSP() const
	{
		return m_allAvailSP;
	}


	/// CPPONLY
	bool allAvailVSP() const
	{
		return m_allAvailVSP;
	}


	vspID resolve(const Population & pop) const;

	/// CPPONLY
	const string & spName() const
	{
		return m_spName;
	}


	/// CPPONLY
	const string & vspName() const
	{
		return m_vspName;
	}


private:
	size_t m_subPop;
	size_t m_virtualSubPop;
	string m_spName;
	string m_vspName;
	bool m_allAvailSP;
	bool m_allAvailVSP;
};

ostream & operator<<(ostream & out, const vspID & vsp);


/** A class to specify (virtual) subpopulation list. Using a dedicated class
 *  allows users to specify a single subpopulation, or a list of (virutal)
 *  subpoulations easily.
 */
class subPopList
{
public:
	// for some unknown reason, std:: is required for this type to be recognized
	// by swig.
	typedef std::vector<vspID> vectorvsp;
	typedef vectorvsp::iterator iterator;
	typedef vectorvsp::const_iterator const_iterator;

public:
	///
	subPopList(PyObject * obj = NULL);

	/// CPPONLY
	subPopList(const vectorvsp & subPops);

	/// CPPONLY
	bool allAvail() const
	{
		return m_allAvail;
	}


	/// CPPONLY
	bool empty() const
	{
		return m_subPops.empty();
	}


	/// CPPONLY
	size_t size() const
	{
		return m_subPops.size();
	}


	size_t __len__() const
	{
		return m_subPops.size();
	}


	/// CPPONLY
	vspID operator[](size_t idx) const
	{
		DBG_FAILIF(idx >= m_subPops.size(), IndexError,
			"Index out of range.");
		return m_subPops[idx];
	}


	/// CPPONLY
	void push_back(const vspID subPop)
	{
		m_subPops.push_back(subPop);
	}


	/// CPPONLY
	bool contains(const vspID subPop) const
	{
		return find(m_subPops.begin(), m_subPops.end(), subPop) != m_subPops.end();
	}


	/// CPPONLY
	bool overlap(const size_t subPop) const
	{
		const_iterator it = m_subPops.begin();
		const_iterator itEnd = m_subPops.end();

		for (; it != itEnd; ++it)
			if (it->subPop() == subPop)
				return true;
		return false;
	}


	/// CPPONLY
	const_iterator begin() const
	{
		return m_subPops.begin();
	}


	/// CPPONLY
	const_iterator end() const
	{
		return m_subPops.end();
	}


	/// CPPONLY
	iterator begin()
	{
		return m_subPops.begin();
	}


	/// CPPONLY
	iterator end()
	{
		return m_subPops.end();
	}


	/// expand ALL_AVAIL and [(ALL_AVAIL, vsp), ...] according to pop
	subPopList expandFrom(const Population & pop) const;

private:
	vectorvsp m_subPops;
	bool m_allAvail;
};


/** This class is the base class of all virtual subpopulation (VSP) splitters,
 *  which provide ways to define groups of individuals in a subpopulation who
 *  share certain properties. A splitter defines a fixed number of named VSPs.
 *  They do not have to add up to the whole subpopulation, nor do they have to
 *  be distinct. After a splitter is assigned to a population, many functions
 *  and operators can be applied to individuals within specified VSPs.
 *
 *  Each VSP has a name. A default name is determined by each splitter but you
 *  can also assign a name to each VSP. The name of a VSP can be retrieved by
 *  function \c BaseVspSplitter.name() or \c Population.subPopName().
 *
 *  Only one VSP splitter can be assigned to a population, which defined VSPs
 *  for all its subpopulations. If different splitters are needed for different
 *  subpopulations, a \c CombinedSplitter can be used.
 */
class BaseVspSplitter
{

public:
	/** This is a virtual class that cannot be instantiated.
	 */
	BaseVspSplitter(const stringList & names = vectorstr()) :
		m_names(names.elems()), m_activated(InvalidValue)
	{
	}


	/** All VSP splitter defines a \c clone() function to create an identical
	 *  copy of itself.
	 */
	virtual BaseVspSplitter * clone() const = 0;

	virtual ~BaseVspSplitter()
	{
	}


	/// Which subpopulation is activated.
	/// CPPONLY
	size_t activatedSubPop() const
	{
		return m_activated;
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	virtual size_t size(const Population & pop, size_t subPop, size_t virtualSubPop) const = 0;

	/** Return the number of VSPs defined by this splitter.
	 */
	virtual size_t numVirtualSubPop() const = 0;

	/** Return \c True if individual \e ind (an index relative to specified
	 *  subpopulation) belongs to specified virtual subpopulation \e vsp.
	 *  CPPONLY
	 */
	virtual bool contains(const Population & pop, size_t ind, vspID vsp) const = 0;

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	virtual void activate(const Population & pop, size_t subPop, size_t virtualSubPop) = 0;

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(size_t subPop)
	{
		(void)subPop;  // avoid a warning message in optimized modules
		DBG_FAILIF(subPop != m_activated, RuntimeError, "Deactivate non-activated virtual subpopulation.");
		m_activated = InvalidValue;
	}


	/** Return the name of VSP \e vsp (an index between \c 0 and
	 *  <tt>numVirtualSubPop()</tt>).
	 */
	virtual string name(size_t vsp) const = 0;

	/** Return the index of a virtual subpopulation from its name. If multiple
	 *  virtual subpopulations share the same name, the first vsp is returned.
	 */
	size_t vspByName(const string & name) const;

protected:
	size_t countVisibleInds(const Population & pop, size_t sp) const;

	vectorstr m_names;

	size_t m_activated;
};

typedef std::vector<BaseVspSplitter *> vectorsplitter;

/** CPPONLY
 *  A functor that judge if an individual belong to a certain virtual subpopulation.
 *  Used by pyIndOperator
 */
class vspFunctor
{
public:
	vspFunctor() : m_pop(NULL), m_splitter(NULL), m_vsp()
	{}

	vspFunctor(const Population & pop, const BaseVspSplitter * splitter, vspID vsp)
		: m_pop(&pop), m_splitter(splitter), m_vsp(vsp)
	{
	}


	bool operator()(size_t ind)
	{
		return m_splitter->contains(*m_pop, ind, m_vsp);
	}


private:
	const Population * m_pop;
	const BaseVspSplitter * m_splitter;
	vspID m_vsp;
};


/** This splitter takes several splitters and stacks their VSPs together. For
 *  example, if the first splitter defines \c 3 VSPs and the second splitter
 *  defines \c 2, the two VSPs from the second splitter become the fourth
 *  (index \c 3) and the fifth (index \c 4) VSPs of the combined splitter.
 *  In addition, a new set of VSPs could be defined as the union of one or more
 *  of the original VSPs. This splitter is usually used to define different
 *  types of VSPs to a population.
 */
class CombinedSplitter : public BaseVspSplitter
{
public:
	/** Create a combined splitter using a list of \e splitters. For example,
	 *  <tt>CombinedSplitter([SexSplitter(), AffectionSplitter()])</tt> defines
	 *  a combined splitter with four VSPs, defined by male (vsp \c 0), female
	 *  (vsp \c 1), unaffected (vsp \c 2) and affected individuals (vsp \c 3).
	 *  Optionally, a new set of VSPs could be defined by parameter \e vspMap.
	 *  Each item in this parameter is a list of VSPs that will be combined to
	 *  a single VSP. For example, <tt>vspMap=[(0, 2), (1, 3)]</tt> in the
	 *  previous example will define two VSPs defined by male or unaffected,
	 *  and female or affected individuals. VSP names are usually determined
	 *  by splitters, but can also be specified using parameter \e names.
	 */
	CombinedSplitter(const vectorsplitter & splitters = vectorsplitter(),
		const intMatrix & vspMap = intMatrix(), const stringList & names = vectorstr());

	/// CPPONLY
	CombinedSplitter(const CombinedSplitter & rhs);

	~CombinedSplitter();

	/// HIDDEN
	BaseVspSplitter * clone() const;

	/// the size of a given virtual subpopulation.
	/// CPPONLY
	size_t size(const Population & pop, size_t subPop, size_t virtualSubPop) const;

	/** Return the number of VSPs defined by this splitter, which is the sum of
	 *  the number of VSPs of all combined splitters.
	 */
	size_t numVirtualSubPop() const
	{
		return m_vspMap.size();
	}


	/** Return \c True if individual \e ind (an index relative to specified
	 *  subpopulation) belongs to specified virtual subpopulation \e vsp.
	 *  CPPONLY
	 */
	bool contains(const Population & pop, size_t ind, vspID vsp) const;

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(const Population & pop, size_t subPop, size_t virtualSubPop);

	/** Return the name of a VSP \e vsp, which is the name a VSP defined by one
	 *  of the combined splitters unless a new set of names is specified. If
	 *  a \e vspMap was used, names from different VSPs will be joined by \c "or".
	 */
	string name(size_t vsp) const;

private:
	/// the splitters
	vector<BaseVspSplitter *> m_splitters;

	/// the splitter correspond to a vsp
	typedef pairu vspPair;
	typedef vector<vspPair> vspList;
	vector<vspList> m_vspMap;
};


/** This splitter takes several splitters and take their intersections as new
 *  VSPs. For example, if the first splitter defines \c 3 VSPs and the second
 *  splitter defines \c 2, \c 6 VSPs will be defined by splitting 3 VSPs
 *  defined by the first splitter each to two VSPs. This splitter is usually
 *  used to define finer VSPs from existing VSPs.
 */
class ProductSplitter : public BaseVspSplitter
{
public:
	/** Create a product splitter using a list of \e splitters. For example,
	 *  <tt>ProductSplitter([SexSplitter(), AffectionSplitter()])</tt> defines
	 *  four VSPs by male unaffected, male affected, female unaffected, and
	 *  female affected individuals. VSP names are usually determined by
	 *  splitters, but can also be specified using parameter \e names.
	 */
	ProductSplitter(const vectorsplitter & splitters = vectorsplitter(),
		const stringList & names = vectorstr());

	/// CPPONLY
	ProductSplitter(const ProductSplitter & rhs);

	~ProductSplitter();

	/// HIDDEN
	BaseVspSplitter * clone() const;

	/// the size of a given virtual subpopulation.
	/// CPPONLY
	size_t size(const Population & pop, size_t subPop, size_t virtualSubPop) const;

	/** Return the number of VSPs defined by this splitter, which is the sum of
	 *  the number of VSPs of all combined splitters.
	 */
	size_t numVirtualSubPop() const
	{
		return m_numVSP;
	}


	/** Return \c True if individual \e ind (an index relative to specified
	 *  subpopulation) belongs to specified virtual subpopulation \e vsp.
	 *  CPPONLY
	 */
	bool contains(const Population & pop, size_t ind, vspID vsp) const;


	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(const Population & pop, size_t subPop, size_t virtualSubPop);

	/** Return the name of a VSP \e vsp, which is the names of indivdual VSPs
	 *  separated by a comma, unless a new set of names is specified for each
	 *  VSP.
	 */
	string name(size_t vsp) const;

private:
	/// the splitters
	vector<BaseVspSplitter *> m_splitters;

	/// total number of vsp
	size_t m_numVSP;

	vector<vectoru> m_subIndexes;
};


/** This splitter defines two VSPs by individual sex. The first VSP consists of
 *  all male individuals and the second VSP consists of all females in a
 *  subpopulation.
 */
class SexSplitter : public BaseVspSplitter
{
public:
	/** Create a sex splitter that defines male and female VSPs. These VSPs
	 *  are named \c Male and \c Female unless a new set of names are specified
	 *  by parameter \e names.
	 */
	SexSplitter(const stringList & names = vectorstr()) : BaseVspSplitter(names)
	{
	}


	/// HIDDEN
	BaseVspSplitter * clone() const
	{
		return new SexSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	size_t size(const Population & pop, size_t subPop, size_t virtualSubPop) const;

	/// Return \c 2.
	size_t numVirtualSubPop() const
	{
		return 2;
	}


	/** Return \c True if individual \e ind (an index relative to specified
	 *  subpopulation) belongs to specified virtual subpopulation \e vsp.
	 *  CPPONLY
	 */
	bool contains(const Population & pop, size_t ind, vspID vsp) const;


	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(const Population & pop, size_t subPop, size_t virtualSubPop);


	/** Return \c "Male" if \e vsp=0 and \c "Female" otherwise, unless a new
	 *  set of names are specified.
	 */
	string name(size_t vsp) const;

};


/** This class defines two VSPs according individual affection status. The
 *  first VSP consists of unaffected invidiauls and the second VSP consists
 *  of affected ones.
 */
class AffectionSplitter : public BaseVspSplitter
{
public:
	/** Create a splitter that defined two VSPs by affection status.These VSPs
	 *  are named \c Unaffected and \c Affected unless a new set of names are
	 *  specified by parameter \e names.
	 */
	AffectionSplitter(const stringList & names = vectorstr()) : BaseVspSplitter(names)
	{
	}


	/// HIDDEN
	BaseVspSplitter * clone() const
	{
		return new AffectionSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	size_t size(const Population & pop, size_t subPop, size_t virtualSubPop) const;

	/// Return 2.
	size_t numVirtualSubPop() const
	{
		return 2;
	}


	/** Return \c True if individual \e ind (an index relative to specified
	 *  subpopulation) belongs to specified virtual subpopulation \e vsp.
	 *  CPPONLY
	 */
	bool contains(const Population & pop, size_t ind, vspID vsp) const;


	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(const Population & pop, size_t subPop, size_t virtualSubPop);


	/** Return \c "Unaffected" if \e vsp=0 and \c "Affected" if \e vsp=1,
	 *  unless a new set of names are specified.
	 */
	string name(size_t vsp) const;

};


/** This splitter defines VSPs according to the value of an information field
 *  of each indivdiual. A VSP is defined either by a value or a range of values.
 */
class InfoSplitter : public BaseVspSplitter
{
public:
	/** Create an infomration splitter using information field \e field. If
	 *  parameter \e values is specified, each item in this list defines a VSP
	 *  in which all individuals have this value at information field \e field.
	 *  If a set of cutoff values are defined in parameter \e cutoff,
	 *  individuals are grouped by intervals defined by these cutoff values.
	 *  For example, <tt>cutoff=[1,2]</tt> defines three VSPs with
	 *  <tt>v < 1</tt>, <tt>1 <= v < 2</tt> and <tt>v >=2</tt> where \c v is the
	 *  value of an individual at information field \e field. If parameter
	 *  \c ranges is specified, each range defines a VSP. For example,
	 *  <tt>ranges=[[1, 3], [2, 5]]</tt> defines two VSPs with
	 *  <tt>1 <= v < 3</tt> and <tt>2 <= 3 < 5</tt>. Of course, only one of the
	 *  parameters \e values, \e cutoff and \e ranges should be defined, and
	 *  values in \e cutoff should be distinct, and in an increasing order. A
	 *  default set of names are given to each VSP unless a new set of names is
	 *  given by parameter \e names.
	 */
	InfoSplitter(string field, const floatList & values = vectorf(),
		const floatList & cutoff = vectorf(), const floatMatrix & ranges = floatMatrix(),
		const stringList & names = vectorstr());

	/// HIDDEN
	BaseVspSplitter * clone() const
	{
		return new InfoSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	size_t size(const Population & pop, size_t subPop, size_t virtualSubPop) const;

	/** Return the number of VSPs defined by this splitter, which is the length
	 *  parameter \e values or the length of \e cutoff plus one, depending on
	 *  which parameter is specified.
	 */
	size_t numVirtualSubPop() const;

	/** Return \c True if individual \e ind (an index relative to specified
	 *  subpopulation) belongs to specified virtual subpopulation \e vsp.
	 *  CPPONLY
	 */
	bool contains(const Population & pop, size_t ind, vspID vsp) const;

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(const Population & pop, size_t subPop, size_t virtualSubPop);


	/** Return the name of a VSP \e vsp, which is <tt>field = value</tt> if VSPs
	 *  are defined by values in parameter \e values, or <tt>field < value</tt>
	 *  (the first VSP), <tt>v1 <= field < v2</tt> and <tt>field >= v</tt> (the
	 *  last VSP) if VSPs are defined by cutoff values. A user-specified name,
	 *  if specified, will be returned instead.
	 */
	string name(size_t vsp) const;

private:
	string m_info;
	//
	vectorf m_values;
	//
	vectorf m_cutoff;
	//
	matrixf m_ranges;
};


/** This splitter divides subpopulations into several VSPs by proportion.
 */
class ProportionSplitter : public BaseVspSplitter
{
public:
	/** Create a splitter that divides subpopulations by \e proportions, which
	 *  should be a list of float numbers (between \c 0 and \c 1) that add up
	 *  to \c 1.  A default set of names are given to each VSP unless a new set
	 *  of names is given by parameter \e names.
	 */
	ProportionSplitter(vectorf const & proportions = vectorf(),
		const stringList & names = vectorstr());

	/// HIDDEN
	BaseVspSplitter * clone() const
	{
		return new ProportionSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	size_t size(const Population & pop, size_t subPop, size_t virtualSubPop) const;

	/** Return the number of VSPs defined by this splitter, which is the length
	 *  of parameter \e proportions.
	 */
	size_t numVirtualSubPop() const;

	/** Return \c True if individual \e ind (an index relative to specified
	 *  subpopulation) belongs to specified virtual subpopulation \e vsp.
	 *  CPPONLY
	 */
	bool contains(const Population & pop, size_t ind, vspID vsp) const;

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(const Population & pop, size_t subPop, size_t virtualSubPop);

	/** Return the name of VSP \e vsp, which is <tt>"Prop p"</tt> where
	 *  <tt>p=propotions[vsp]</tt>. A user specified name will be returned if
	 *  specified.
	 */
	string name(size_t vsp) const;

private:
	vectorf m_proportions;
};


/** This class defines a splitter that groups individuals in certain ranges
 *  into VSPs.
 */
class RangeSplitter : public BaseVspSplitter
{
public:
	/** Create a splitter according to a number of individual ranges defined
	 *  in \e ranges. For example,
	 *  <tt>RangeSplitter(ranges=[[0, 20], [40, 50]])</tt> defines two VSPs.
	 *  The first VSP consists of individuals \c 0, \c 1, ..., \c 19, and the
	 *  sceond VSP consists of individuals \c 40, \c 41, ..., \c 49. Note that
	 *  a nested list has to be used even if only one range is defined. A
	 *  default set of names are given to each VSP unless a new set of names is
	 *  given by parameter \e names.
	 */
	RangeSplitter(const intMatrix & ranges, const stringList & names = vectorstr());

	/// HIDDEN
	BaseVspSplitter * clone() const
	{
		return new RangeSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	size_t size(const Population & pop, size_t subPop, size_t virtualSubPop) const;

	/** Return the number of VSPs, which is the number of ranges defined in
	 *  parameter \e ranges.
	 */
	size_t numVirtualSubPop() const;

	/** Return \c True if individual \e ind (an index relative to specified
	 *  subpopulation) belongs to specified virtual subpopulation \e vsp.
	 *  CPPONLY
	 */
	bool contains(const Population & pop, size_t ind, vspID vsp) const;

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(const Population & pop, size_t subPop, size_t virtualSubPop);

	/** Return the name of VSP \e vsp, which is <tt>"Range [a, b)"</tt> where
	 *  <tt>[a, b)</tt> is range <tt>ranges[vsp]</tt>. A user specified name
	 *  will be returned if specified.
	 */
	string name(size_t vsp) const;

private:
	matrixi m_ranges;
};


/** This class defines a VSP splitter that defines VSPs according to individual
 *  genotype at specified loci.
 */
class GenotypeSplitter : public BaseVspSplitter
{
public:
	/** Create a splitter that defines VSPs by individual genotype at \e loci
	 *  (can be indexes or names of one or more loci). Each list in a list
	 *  \e allele defines a VSP, which is a list of allowed alleles at these
	 *  \e loci. If only one VSP is defined, the outer list of the nested list
	 *  can be ignored. If phase if true, the order of alleles in each list is
	 *  significant. If more than one set of alleles are given, Individuals
	 *  having either of them is qualified.
	 *
	 *  For example, in a haploid population, <tt>loci=1, alleles=[0, 1]</tt>
	 *  defines a VSP with individuals having allele \c 0 or \c 1 at locus \c 1,
	 *  <tt>alleles=[[0, 1], [2]]</tt> defines two VSPs with indivdiuals in the
	 *  second VSP having allele \c 2 at locus \c 1. If multiple loci are
	 *  involved, alleles at each locus need to be defined. For example,
	 *  VSP defined by <tt>loci=[0, 1], alleles=[0, 1, 1, 1]</tt> consists of
	 *  individuals having alleles <tt>[0, 1]</tt> or <tt>[1, 1]</tt> at loci
	 *  <tt>[0, 1]</tt>.
	 *
	 *  In a haploid population, <tt>loci=1, alleles=[0, 1]</tt> defines a VSP
	 *  with individuals having genotype <tt>[0, 1]</tt> or <tt>[1, 0]</tt> at
	 *  locus \c 1. <tt>alleles[[0, 1], [2, 2]]</tt> defines two VSPs with
	 *  indivdiuals in the second VSP having genotype <tt>[2, 2]</tt> at locus
	 *  \c 1. If \e phase is set to \c True, the first VSP will only has
	 *  individuals with genotype <tt>[0, 1]</tt>. In the multiple loci case,
	 *  alleles should be arranged by haplotypes, for example,
	 *  <tt>loci=[0, 1], alleles=[0, 0, 1, 1], phase=True</tt> defines a VSP with
	 *  individuals having genotype <tt>-0-0-, -1-1-</tt> at loci \c 0 and \c 1.
	 *  If <tt>phase=False</tt> (default), genotypes <tt>-1-1-, -0-0-</tt>,
	 *  <tt>-0-1-</tt> and <tt>-1-0-</tt> are all allowed.
	 *
	 *  A default set of names are given to each VSP unless a new set of names
	 *  is given by parameter \e names.
	 */
	GenotypeSplitter(const lociList & loci,
		const intMatrix & alleles, bool phase = false,
		const stringList & names = vectorstr());

	/// HIDDEN
	BaseVspSplitter * clone() const
	{
		return new GenotypeSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	size_t size(const Population & pop, size_t subPop, size_t virtualSubPop) const;

	/// number of virtual subpops of subpopulation sp
	size_t numVirtualSubPop() const;

	/** Return \c True if individual \e ind (an index relative to specified
	 *  subpopulation) belongs to specified virtual subpopulation \e vsp.
	 *  CPPONLY
	 */
	bool contains(const Population & pop, size_t ind, vspID vsp) const;

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(const Population & pop, size_t subPop, size_t virtualSubPop);

	/** Return name of VSP \e vsp, which is <tt>"Genotype loc1,loc2:genotype"</tt>
	 *  as defined by parameters \e loci and \e alleles. A user provided name
	 *  will be returned if specified.
	 */
	string name(size_t vsp) const;

private:
	bool match(const Individual * ind, const vectori & alleles) const;

	bool matchSingle(const Individual * ind, const vectori & alleles) const;

private:
	lociList m_loci;
	matrixi m_alleles;
	bool m_phase;
};

}
#endif
