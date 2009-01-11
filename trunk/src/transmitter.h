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

#ifndef _TRANSMITTER_H
#define _TRANSMITTER_H
/**
   \file
   \brief head file of class recombinator:public baseOperator
 */
#include "operator.h"

#include <iterator>
using std::ostream;
using std::ostream_iterator;

namespace simuPOP {

/** This during mating operator is the base class of all genotype transmitters.
 *  It is made available to users because it provides a few member functions
 *  that can be used by derived transmitters, and by customized Python
 *  during mating operators.
 */
class genoTransmitter : public baseOperator
{
public:
	/** Create a base genotype transmitter.
	 */
	genoTransmitter(int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(),
		const vectorstr & infoFields = vectorstr()) :
		baseOperator("", DuringMating, begin, end, step, at, rep, subPop, infoFields),
		m_ploidy(0), m_hasCustomizedChroms(false), m_lociToCopy(0), m_chromIdx(0)
	{
		setFormOffGenotype(true);
	}


	/// Deep copy of a base genotype transmitter.
	baseOperator * clone() const
	{
		return new genoTransmitter(*this);
	}


	/** Clear (set alleles to zero) chromosome \e chrom on the \e ploidy-th
	 *  homologous set of chromosomes of individual \e ind.
	 */
	void clearChromosome(const individual & ind, int ploidy, int chrom);

	/** Transmit chromosome \e chrom on the \e parPloidy set of homologous
	 *  chromosomes from \e parent to the \e ploidy set of homologous
	 *  chromosomes of \e offspring.
	 */
	void copyChromosome(const individual & parent, int parPloidy,
		individual & offspring, int ploidy, int chrom);

	/** Transmit the \e parPloidy set of homologous chromosomes from \e parent
	 *  to the \e ploidy set of homologous chromosomes of \e offspring.
	 *  Customized chromosomes are not copied.
	 */
	void copyChromosomes(const individual & parent, int parPloidy,
		individual & offspring, int ploidy);

	virtual string __repr__()
	{
		return "<simuPOP::genoTransmitter>" ;
	}

	/** Initialize a base genotype operator for a population. This function should be
	 *  called before any other functions are used to transmit genotype.
	 */
	void initialize(const population & pop);

	/// CPPONLY
	bool applyDuringMating(population & pop,
	                       RawIndIterator offspring,
	                       individual * dad = NULL,
	                       individual * mom = NULL)
	{
		throw SystemError("The base genotype transmitter does not provide any function to transmit genotype");
	}


protected:
	// cache some genostructor information for
	// faster performance
	UINT m_ploidy;
	bool m_hasCustomizedChroms;
	vectoru m_lociToCopy;
	vectoru m_chromIdx;
};

/** This during mating operator copies parental genotype directly to offspring.
 *  This operator works for all mating schemes when one or two parents are
 *  involved. If both parents are passed, maternal genotype are copied.
 */
class cloneGenoTransmitter : public genoTransmitter
{
public:
	/** Create a clone genotype transmitter.
	 */
	cloneGenoTransmitter(int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(),
		const vectorstr & infoFields = vectorstr()) :
		genoTransmitter(begin, end, step, at, rep, subPop, infoFields)
	{
		setFormOffGenotype(true);
	}

	/// Deep copy of a clone genotype transmitter.
	baseOperator * clone() const
	{
		return new cloneGenoTransmitter(*this);
	}


	/// HIDDEN
	virtual string __repr__()
	{
		return "<simuPOP::cloneGenoTransmitter>" ;
	}


	/// CPPONLY
	bool applyDuringMating(population & pop,
		RawIndIterator offspring,
		individual * dad = NULL,
		individual * mom = NULL);

};


/** Mendelian offspring generator accepts two parents and pass their
   genotype to a number of offspring following Mendelian's law. Basically,
   one of the paternal chromosomes is chosen randomly to form the paternal
   copy of the offspring, and one of the maternal chromosome is chosen
   randomly to form the maternal copy of the offspring. The number of offspring
   produced is controled by parameters \c numOffspring, \c numOffspringFunc,
   \c maxNumOffspring and \c mode. Recombination will not happen unless
   a during-mating operator recombinator is used.
 */
class mendelianGenoTransmitter : public genoTransmitter
{
public:
	/** Create a Mendelian genotype transmitter.
	 */
	mendelianGenoTransmitter(int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(),
		const vectorstr & infoFields = vectorstr()) :
		genoTransmitter(begin, end, step, at, rep, subPop, infoFields),
		m_chromX(-1), m_chromY(-1), m_numChrom(0)
	{
	}

	/// Deep copy of a Mendelian genotype transmitter.
	baseOperator * clone() const
	{
		return new mendelianGenoTransmitter(*this);
	}


	/// HIDDEN
	virtual string __repr__()
	{
		return "<simuPOP::mendelianGenoTransmitter>" ;
	}

	/// CPPONLY
	virtual bool applyDuringMating(population & pop,
		RawIndIterator offspring,
		individual * dad = NULL,
		individual * mom = NULL);

	/** Initialize a base genotype operator for a population. This function should be
	 *  called before function \c transmitGenotype is used to transmit genotype.
	 */
	void initialize(const population & pop);

	/** Transmit genotype from parent to offspring, and fill the \e ploidy
	 *  homologous set of chromosomes. This function does not set genotypes of
	 *  customized chromosomes and handles sex chromosomes properly, according
	 *  to offspring sex and \c ploidy.
	 */
	void transmitGenotype(const individual & parent,
		individual & offspring, int ploidy);

protected:
	// cache chromBegin, chromEnd for better performance.
	vectoru m_chIdx;
	int m_chromX;
	int m_chromY;
	UINT m_numChrom;
};


/** selfing offspring generator works similarly as a mendelian offspring
   generator but a single parent produces both the paternal and maternal
   copy of the offspring chromosomes. This offspring generator accepts a
   dipload parent. A random copy of the parental chromosomes is chosen
   randomly to form the parental copy of the offspring chromosome, and
   is chosen randomly again to form the maternal copy of the offspring
   chromosome.
 */
class selfingGenoTransmitter : public mendelianGenoTransmitter
{
public:
	/// Create a self-fertilization genotype transmitter.
	selfingGenoTransmitter(int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(),
		const vectorstr & infoFields = vectorstr())
		: mendelianGenoTransmitter(begin, end, step, at, rep, subPop, infoFields)
	{
	}


	/// Deep copy of a selfing genotype transmitter
	baseOperator * clone() const
	{
		return new selfingGenoTransmitter(*this);
	}


	/// HIDDEN
	virtual string __repr__()
	{
		return "<simuPOP::selfingGenoTransmitter>" ;
	}

	/// CPPONLY
	bool applyDuringMating(population & pop,
		RawIndIterator offspring,
		individual * dad = NULL,
		individual * mom = NULL);

};


/** haplodiploid offspring generator mimics sex-determination in honey bees.
   Given a female (queen) parent and a male parent, the female is considered
   as diploid with two set of chromosomes, and the male is condiered as haploid.
   Actually, the first set of male chromosomes are used. During mating,
   female produce eggs, subject to potential recombination and gene conversion,
   while male sperm is identical to the parental chromosome.

   Female offspring has two sets of chromosomes, one from mother and one from
   father. Male offspring has one set of chromosomes from his mother.
   <applicability>haplodiploid only</applicability>
 */
class haplodiploidGenoTransmitter : public mendelianGenoTransmitter
{
public:
	/// Create a haplodiploid genotype transmitter.
	haplodiploidGenoTransmitter(int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(),
		const vectorstr & infoFields = vectorstr())
		: mendelianGenoTransmitter(begin, end, step, at, rep, subPop, infoFields),
		m_copier()
	{
	}


	/// Deep copy of a haplodiploid transmitter.
	baseOperator * clone() const
	{
		return new haplodiploidGenoTransmitter(*this);
	}

	/// HIDDEN
	virtual string __repr__()
	{
		return "<simuPOP::haplodiploidGenoTransmitter>" ;
	}

	/// HIDDEN
	void initialize(const population & pop);

	/// CPPONLY
	virtual bool applyDuringMating(population & pop,
		RawIndIterator offspring,
		individual * dad = NULL,
		individual * mom = NULL);

private:
	genoTransmitter m_copier;
};


/** This geno transmitter assumes that the first homologous copy of several (or
 *  all) Customized chromosomes are copies of mitochondrial chromosomes. It
 *  transmits these chromosomes randomly from the female parent.
 */
class mitochondrialGenoTransmitter : public genoTransmitter
{
public:
	/** Createa a mitochondrial genotype transmitter that treats all Customized
	 *  chromosomes, or a list of chromosomes specified by \e chroms, as human
	 *  mitochondrial chromosomes. It transmits these chromosomes randomly from
	 *  the female parent to offspring of both sexes.
	 */
	mitochondrialGenoTransmitter(const vectoru & chroms = vectoru(),
		int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(),
		const vectorstr & infoFields = vectorstr())
		: genoTransmitter(begin, end, step, at, rep, subPop, infoFields),
		m_chroms(chroms), m_mitoChroms(0), m_numLoci(0)
	{
		// this is intended to be an auxillary genotype transmitter.
		setFormOffGenotype(false);
	}

	/// Deep copy of a mitochondrial genotype transmitter.
	baseOperator * clone() const
	{
		return new mitochondrialGenoTransmitter(*this);
	}

	/// HIDDEN
	virtual string __repr__()
	{
		return "<simuPOP::mitochondrialGenoTransmitter>" ;
	}

	/// HIDDEN
	void initialize(const population & pop);

	/// CPPONLY
	virtual bool applyDuringMating(population & pop,
		RawIndIterator offspring,
		individual * dad = NULL,
		individual * mom = NULL);

private:
	// this is user input.
	vectoru m_chroms;

	// this is the temporary holder for different populaitons
	vectoru m_mitoChroms;

	//
	UINT m_numLoci;
};


/**
   In simuPOP, only one recombinator is provided. Recombination events between loci
   a/b and b/c are independent, otherwise there will be some linkage between loci. Users
   need to specify physical recombination rate between adjacent loci. In addition,
   for the recombinator
   \li it only works for diploid (and for females in haplodiploid) populations.
   \li the recombination rate must be comprised between \c 0.0 and \c 0.5. A recombination
   rate of \c 0.0 means that the loci are completely linked, and thus behave together
   as a single linked locus. A recombination rate of \c 0.5 is equivalent to free of
   recombination. All other values between \c 0.0 and \c 0.5 will represent various
   linkage intensities	between adjacent pairs of loci. The recombination rate is
   equivalent to <tt>1-linkage</tt> and represents the probability that the allele
   at the next locus is randomly drawn.
   \li it works for selfing. I.e., when only one parent is provided, it will be
   recombined twice, producing both maternal and paternal chromosomes of the
   offspring.
   \li conversion is allowed. Note that conversion will nullify many recombination
    events, depending on the parameters chosen.
 */

class recombinator : public genoTransmitter
{
public:
	/// recombine chromosomes from parents
	/**
	   \param intensity intensity of recombination. The actual recombination rate
	    between two loci is determined by <tt>intensity*locus distance (between them)</tt>.
	   \param rate recombination rate regardless of locus distance after all \c afterLoci.
	    It can also be an array of recombination rates. Should have the same length
	    as \c afterLoci or \c totNumOfLoci().
	    The recombination rates are independent of locus distance.
	   \param afterLoci an array of locus indexes. Recombination will occur after these
	    loci. If \c rate is also specified, they should have the same length. Default
	    to all loci (but meaningless for those loci located at the end of a chromosome).
	    If this parameter is given, it should be ordered, and can not include loci at
	   the end of a chromosome.

	 * \e convMode can take the following forms
	 * NoConversion: no conversion
	 * (NumMarkers, prob, n): converts a fixed number of markers
	 * (GeometricDistribution, prob, p): An geometric distribution is used to
	        determine how many markers will be converted.
	 * (TractLength, prob, n): converts a fixed length of tract.
	 * (ExponentialDistribution, prob, p): An exponential distribution with parameter
	        \c convLen will be used to determine track length.
	 * The first number is that probability of conversion event among all recombination
	    events. When a recombination event happens, it may become a recombination event
	    if the Holliday junction is resolved/repaired successfully, or a
	    conversion event if the junction is not resolved/repaired. The
	    default \c convProb is 0, meaning no conversion event at all.
	    Note that the ratio of conversion to recombination events varies greatly from
	    study to study, ranging from 0.1 to 15 (Chen et al, Nature Review Genetics, 2007).
	    This translate to 0.1/0.9~0.1 to 15/16~0.94 of this parameter. When

	    Note that
	   \li conversion tract length is usually short, and is estimated to be
	        between 337 and 456 bp, with overall range between maybe 50 - 2500 bp.
	   \li simuPOP does not impose a unit for marker distance so your choice
	        of \c convParam needs to be consistent with your unit. In the HapMap dataset,
	        cM is usually assumed and marker distances are around 10kb (0.001cM ~- 1kb).
	        Gene conversion can largely be ignored. This is important when
	        you use distance based conversion mode such as \c CONVERT_TrackLength or
	   \c CONVERT_ExponentialDistribution.
	   \li After a track length is determined, if a second recombination
	        event happens within this region, the track length will be shortened.
	        Note that conversion is identical to double recombination under
	        this context.
	   \param haplodiploid If set to true, the first copy of paternal chromosomes
	        is copied directly as the paternal chromosomes of the offspring. This
	        is because haplodiploid male has only one set of chromosome.

	   \note There is no recombination between sex chromosomes of male individuals
	   if <tt>sexChrom()=True</tt>. This may change later if the exchanges
	   of genes between pseudoautosomal regions of \c XY need to be modeled.

	 */
	recombinator(double intensity = -1, vectorf rate = vectorf(), vectoru loci = vectoru(),
		const floatList & convMode = NoConversion,
		int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(),
		const vectorstr & infoFields = vectorstr());


	/// deep copy of a recombinator
	virtual baseOperator * clone() const
	{
		return new recombinator(*this);
	}


	/// used by Python print function to print out the general information of the recombinator
	virtual string __repr__()
	{
		return "<simuPOP::recombination>" ;
	}


	/// return recombination counts (only valid in standard modules)
	UINT recCount(size_t idx)
	{
		DBG_FAILIF(idx >= m_recCount.size(), IndexError,
			"RecCount index out of range");
#ifndef OPTIMIZED
		return m_recCount[idx];
#else
		return 0;
#endif
	}


	/// return recombination counts (only valid in standard modules)
	vectoru recCounts()
	{
#ifndef OPTIMIZED
		return m_recCount;
#else
		return vectoru();
#endif
	}


	/// return the count of conversion of a certain size (only valid in standard modules)
	ULONG convCount(size_t size)
	{
#ifndef OPTIMIZED
		return m_convSize[size];
#else
		return 0;
#endif
	}


	/// return the count of conversions of all sizes (only valid in standard modules)
	std::map<int, int> convCounts()
	{
#ifndef OPTIMIZED
		return m_convSize;
#else
		return std::map<int, int>();
#endif
	}


	void initialize(const population & pop);

	// this function implement how to recombine
	// parental chromosomes and set one copy of offspring chromosome
	// bt contains the bernulli trailer

	void transmitGenotype(const individual & parent,
		individual & offspring, int ploidy);

	/// apply the recombinator during mating
	/// CPPONLY
	virtual bool applyDuringMating(population & pop,
		RawIndIterator offspring,
		individual * dad, individual * mom);

private:
	/// determine number of markers to convert
	int markersConverted(size_t index, const individual & ind);

private:
	/// intensity
	double m_intensity;

	/// differnt rates
	vectorf m_rate;

	/// initial parameter
	vectoru m_afterLoci;

	/// position to recombine, changed to fit a special pop
	vectoru m_recBeforeLoci;

	floatList m_convMode;

	/// bernulli trials
	//  vector<BernulliTrials*> m_bt;
	BernulliTrials m_bt;

	// locataion of special chromosomes
	int m_chromX;
	int m_chromY;
	int m_customizedBegin;
	int m_customizedEnd;

#ifndef OPTIMIZED
	/// report the number of recombination events
	vectoru m_recCount;

	// report the tract length of conversions
	std::map<int, int> m_convSize;
#endif
	/// algorithm to use (frequent or seldom recombinations)
	int m_algorithm;
};

}
#endif
