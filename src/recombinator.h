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
/// recombination and conversion
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
			of convLen needs to be consistent with your unit. In the HapMap dataset,
			cM is usually assumed and marker distances are around 10kb (0.001cM ~- 1kb).
			Gene conversion can largely be ignored. This is important when
			you use distance based conversion mode such as \c CONVERT_TrackLength or
			\c CONVERT_ExponentialDistribution.
		\li After a track length is determined, if a second recombination
			event happens within this region, the track length will be shortened.
			Note that conversion is identical to double recombination under
			this context.

	 \note There is no recombination between sex chromosomes of male individuals
	   if <tt>sexChrom()=True</tt>. This may change later if the exchanges
	   of genes between pseudoautosomal regions of \c XY need to be modeled.

	 \test src_recombinator.log Operator \c recombinator
	 */
	recombinator(double intensity = -1,
	             vectorf rate = vectorf(),
	             vectoru afterLoci = vectoru(),
	             double maleIntensity = -1,
	             vectorf maleRate = vectorf(),
	             vectoru maleAfterLoci = vectoru(),
				 double convProb = 0, // no conversion
				 UINT convMode = CONVERT_NumMarkers,
				 double convParam = 1.,
	             int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
	             int rep = REP_ALL, int grp = GRP_ALL, const vectorstr & infoFields = vectorstr())
		:
		baseOperator("", "", DuringMating, begin, end, step, at, rep, grp, infoFields)
		, m_intensity(intensity), m_maleIntensity(maleIntensity),
		m_rate(rate), m_maleRate(maleRate),
		m_afterLoci(afterLoci), m_maleAfterLoci(maleAfterLoci),
		m_recBeforeLoci(0), m_maleRecBeforeLoci(0),
		m_convProb(convProb), m_convMode(convMode), m_convParam(convParam), 
		m_bt(rng()), m_maleBt(rng()), m_recCount(0), m_algorithm(0)
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


	/// return recombination count
	ULONG recCount(size_t locus)
	{
		DBG_ASSERT(locus < m_recCount.size(), IndexError,
			"locus index " + toStr(locus) + " is out of range");
		return m_recCount[locus];
	}


	/// return recombination counts
	vectoru recCounts()
	{
		return m_recCount;
	}


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
		individual * parent,                                                        // one of the parent
		RawIndIterator & offspring,                                                 // offspring
		int offPloidy,                                                              // which offspring ploidy to fill
		BernulliTrials & bt,
		const vectoru & recBeforeLoci,
		bool setSex = false);
	
	/// determine number of markers to convert
	int markersConverted(size_t index, individual * ind);

	/// this function takes intensity, rate, afterLoci, ...
	/// inputs and return a bernulli trailer and a recBeforeLoci
	/// vector.
	void prepareRecRates(population & pop,
		double intensity,
		vectorf rate,
		vectoru afterLoci,                                                  //
		bool sexChrom,                                                      // whether or not recombine the last chromosome
		vectoru & recBeforeLoci,                                            // return before loci vector
		vectorf & vecP);                                                    // return recombination rate

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

	/// whether or not set sex (population having sex chromosome)
	bool m_hasSexChrom;

	/// report the number of recombination events
	vectoru m_recCount;

	/// algorithm to use (frequent or seldom recombinations)
	int m_algorithm;
};

}
#endif
