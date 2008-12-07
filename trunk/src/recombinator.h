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

#ifndef _RECOMBINATOR_H
#define _RECOMBINATOR_H
/**
   \file
   \brief head file of class recombinator:public baseOperator
 */
#include "operator.h"

#include <iterator>
using std::ostream;
using std::ostream_iterator;

namespace simuPOP {
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

class recombinator : public baseOperator
{
public:
#define CONVERT_NumMarkers                   1
#define CONVERT_TractLength                  2
#define CONVERT_ExponentialDistribution      3
#define CONVERT_GeometricDistribution        4

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
	   \param maleIntensity recombination intensity for male individuals. If given,
	   parameter \c intensity will be considered as female intensity.
	   \param maleRate recombination rate for male individuals. If given,
	   parameter \c rate will be considered as female recombination rate.
	   \param maleAfterLoci if given, males will recombine at different locations.
	   \param convProb The probability of conversion event among all recombination
	    events. When a recombination event happens, it may become a recombination event
	    if the Holliday junction is resolved/repaired successfully, or a
	    conversion event if the junction is not resolved/repaired. The
	    default \c convProb is 0, meaning no conversion event at all.
	    Note that the ratio of conversion to recombination events varies greatly from
	    study to study, ranging from 0.1 to 15 (Chen et al, Nature Review Genetics, 2007).
	    This translate to 0.1/0.9~0.1 to 15/16~0.94 of this parameter. When
	   \c convProb is 1, all recombination events will be conversion events.
	   \param convMode conversion mode, determines how track length is determined.
	   \li CONVERT_NumMarkers Converts a fixed number of markers.
	   \li CONVERT_GeometricDistribution An geometric distribution is used to
	        determine how many markers will be converted.
	   \li CONVERT_TractLength Converts a fixed length of tract.
	   \li CONVERT_ExponentialDistribution An exponential distribution with parameter
	   \c convLen will be used to determine track length.
	   \param convParam Parameter for the conversion process. The exact meaning of this
	    parameter is determined by \c convMode. Note that
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
	recombinator(double intensity = -1,
		vectorf rate = vectorf(),
		vectoru afterLoci = vectoru(),
		double maleIntensity = -1,
		vectorf maleRate = vectorf(),
		vectoru maleAfterLoci = vectoru(),
		double convProb = 0,          // no conversion
		UINT convMode = CONVERT_NumMarkers,
		double convParam = 1.,
		int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr())
		:
		baseOperator("", "", DuringMating, begin, end, step, at, rep, subPop, infoFields)
		, m_intensity(intensity), m_maleIntensity(maleIntensity),
		m_rate(rate), m_maleRate(maleRate),
		m_afterLoci(afterLoci), m_maleAfterLoci(maleAfterLoci),
		m_recBeforeLoci(0), m_maleRecBeforeLoci(0),
		m_convProb(convProb), m_convMode(convMode), m_convParam(convParam),
		m_bt(rng()), m_maleBt(rng()),
#ifndef OPTIMIZED
		m_recCount(0), m_convSize(),
#endif
		m_algorithm(0)
	{
		// tells mating schemes that this operator will form
		// the genotype of offspring so they do not have to
		// generate default genotype for offspring
		this->setFormOffGenotype(true);

		DBG_FAILIF(fcmp_lt(m_convProb, 0) || fcmp_gt(m_convProb, 1),
			ValueError, "Conversion probability should be between 0 and 1");
	};

	virtual ~recombinator()
	{
	}


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


	/// return recombination count at a locus (only valid in standard modules)
	ULONG recCount(size_t locus)
	{
#ifndef OPTIMIZED
		DBG_ASSERT(locus < m_recCount.size(), IndexError,
			"locus index " + toStr(locus) + " is out of range");
		return m_recCount[locus];
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

	/** Recombine parental chromosomes of \e parent	and pass them to offspring
	 *  \e off. The homologous chromosomes ofe parent will be recombined
	 *  twice and form both homologous sets of the offspring, as if \e parent
	 *  mates with itself (a selfing inheritance model). If sex chromosomes
	 *  are present, offspring sex will be determined by which sex chromosomes
	 *  are inherited by \e off. Random sex is assigned to \e off otherwise.
	 */
	void produceOffspring(const individual & parent, individual & off);

	/** Recombine parental chromosomes and pass them to offspring \e off. A
	 *  Mendelian inheritance model will be used, which recombine homologous
	 *  sets of chromosomes of \e mom and \e dad and pass them as the first and
	 *  second sets of homologous chromosomes to offspring \e off, respectively.
	 *  If sex chromosomes are present, offspring sex is determined by which
	 *  sex chromosomes are inherited by \e off. Random sex is assigned to
	 *  \e off otherwise.
	 */
	void produceOffspring(const individual & mom, const individual & dad,
		individual & off);

	/// apply the recombinator during mating
	/// CPPONLY
	virtual bool applyDuringMating(population & pop,
		RawIndIterator offspring,
		individual * dad = NULL,
		individual * mom = NULL);

private:
	// this function implement how to recombine
	// parental chromosomes and set one copy of offspring chromosome
	// bt contains the bernulli trailer
	void recombine(
		const individual & parent,                                                          // one of the parent
		individual & offspring,                                                             // offspring
		int offPloidy,                                                                      // which offspring ploidy to fill
		BernulliTrials & bt,
		const vectoru & recBeforeLoci,
		bool setSex = false);

	/// determine number of markers to convert
	int markersConverted(size_t index, const individual & ind);

	/// this function takes intensity, rate, afterLoci, ...
	/// inputs and return a bernulli trailer and a recBeforeLoci
	/// vector.
	void prepareRecRates(const population & pop,
		double intensity,
		vectorf rate,
		vectoru afterLoci,                                                  //
		vectoru & recBeforeLoci,                                            // return before loci vector
		vectorf & vecP,
		Sex sex);                                                           // return recombination rate

	void copyParentalGenotype(const individual & parent,
		individual & it, int ploidy);

private:
	/// intensity
	double m_intensity;
	double m_maleIntensity;

	/// differnt rates
	vectorf m_rate;
	vectorf m_maleRate;

	/// initial parameter
	vectoru m_afterLoci;
	vectoru m_maleAfterLoci;

	/// position to recombine, changed to fit a special pop
	vectoru m_recBeforeLoci;
	vectoru m_maleRecBeforeLoci;

	double m_convProb;

	UINT m_convMode;

	double m_convParam;

	/// bernulli trials
	//  vector<BernulliTrials*> m_bt;
	BernulliTrials m_bt, m_maleBt;

	// locataion of special chromosomes
	int m_chromX;
	int m_chromY;
	vectoru m_customized;

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
