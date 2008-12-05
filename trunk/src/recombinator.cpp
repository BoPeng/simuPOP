/***************************************************************************
*   Copyright (C) 2004 by Bo Peng										 *
*   bpeng@rice.edu														*
*																		 *
*   $LastChangedDate$
*   $Rev$
*
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or	 *
*   (at your option) any later version.								   *
*																		 *
*   This program is distributed in the hope that it will be useful,	   *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of		*
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the		 *
*   GNU General Public License for more details.						  *
*																		 *
*   You should have received a copy of the GNU General Public License	 *
*   along with this program; if not, write to the						 *
*   Free Software Foundation, Inc.,									   *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.			 *
***************************************************************************/

#include "recombinator.h"

using std::min;
using std::max;

namespace simuPOP {
void recombinator::prepareRecRates(const population & pop,
                                   double intensity,
                                   vectorf rate,
                                   vectoru afterLoci,                   //
                                   vectoru & recBeforeLoci,             // return before loci vector
                                   vectorf & vecP,                      // return recombination rate
                                   Sex sex)
{
	if (!recBeforeLoci.empty() )
		return;

	DBG_FAILIF(intensity < 0 && rate.empty(), ValueError,
		"You should specify intensity, or rate "
		"(a number or a sequence of recombination rates.)");

	DBG_FAILIF(rate.size() > 1 && afterLoci.empty(), ValueError,
		"When more than one rates are given, afterLoci should be"
		" explicitly specified.");

	DBG_FAILIF(rate.size() > 1 && rate.size() != afterLoci.size(),
		ValueError, "If both rates and atLoci are specified, "
		            "they should have the same length.");

	bool useLociDist = rate.empty();

	recBeforeLoci.clear();
	vecP.clear();
	for (UINT ch = 0; ch < pop.numChrom(); ++ch) {
		UINT chBegin = pop.chromBegin(ch);
		UINT chEnd = pop.chromEnd(ch);

		if (std::find(m_mitochondrial.begin(), m_mitochondrial.end(), ch) != m_mitochondrial.end() ||  // mitochondrial does not recombine
		    (sex == Male && (static_cast<int>(ch) == m_chromX || static_cast<int>(ch) == m_chromY)) ||
		    (sex == Female && static_cast<int>(ch) == m_chromY)) {
			DBG_DO(DBG_RECOMBINATOR, cout << "Ignoring chromosome " << ch
				<< " Mitochondrial: " << m_mitochondrial << " Chrom X: " << m_chromX
				<< " Chrom Y: " << m_chromY << endl);
			continue;
		}

		if (afterLoci.empty()) {
			// get loci distance * rate and then recombinant points
			for (UINT loc = chBegin; loc < chEnd - 1; ++loc) {
				recBeforeLoci.push_back(loc + 1);
				double r = useLociDist ? ((pop.locusPos(loc + 1) - pop.locusPos(loc)) * intensity) : rate[0];

				DBG_WARNING(fcmp_gt(r, 0.5),
					"Recombination rate after marker " + toStr(loc) + " is out of range ("
					+ toStr(r) + " ) so it is set to 0.5. This may happen \n"
					             "when you use recombination intensity instead of rate, and your loci \n"
					             "distance is too high.)");
				vecP.push_back(min(0.5, r));
			}
		} else {                                                                          // afterLoci not empty
			DBG_FAILIF(rate.size() > 1 && rate.size() != afterLoci.size(), SystemError,
				"If an array is given, rates and afterLoci should have the same length");

			// get loci distance * rate and then recombinant points
			for (UINT loc = chBegin; loc < chEnd - 1; ++loc) {
				// if this locus will be recombined.
				vectoru::iterator pos = find(afterLoci.begin(), afterLoci.end(), loc);
				if (pos != afterLoci.end()) {
					double r = 0;
					if (useLociDist)
						r = intensity > 0 ? ((pop.locusPos(loc + 1) - pop.locusPos(loc)) * intensity) : r;
					else if (rate.size() == 1 && !useLociDist)
						r = max(rate[0], 0.);
					else
						r = rate[pos - afterLoci.begin()];
					recBeforeLoci.push_back(loc + 1);
					vecP.push_back(r);

					DBG_ASSERT(fcmp_ge(vecP[vecP.size() - 1], 0) && fcmp_le(vecP[vecP.size() - 1], 1),
						ValueError,
						"Recombination rate should be in [0,1]. (Maybe your loci distance is too high.)");
				}
			}
		}
		// after each chromosome ...
		recBeforeLoci.push_back(chEnd);
		vecP.push_back(0.5);
	}
	DBG_DO(DBG_RECOMBINATOR, cout << "Specify after Loci. With rates "
								  << vecP << " before " << recBeforeLoci << endl);

	DBG_ASSERT(vecP.size() == recBeforeLoci.size(), SystemError,
		"Rate and before loci should have the same length.");

	DBG_ASSERT(recBeforeLoci.back() == pop.totNumLoci(),
		SystemError, "The last beforeLoci elem should be total number of loci.");

	DBG_ASSERT(vecP.back() == .5, SystemError,
		"The last elem of rate should be half.");

	// initialize recombination counter,
	// This will count recombination events after
	// each locus.
	DBG_DO_(m_recCount.resize(pop.totNumLoci(), 0));
	return;
}


int recombinator::markersConverted(size_t index, const individual & ind)
{
	// IMPORTANT: if conversion length reaches end of chromosome
	// this is an recombination! Otherwise, conversion will
	// interfere with free crossover between chromosomes
	if (m_convMode == CONVERT_NumMarkers || m_convMode == CONVERT_GeometricDistribution) {
		UINT num = 0;
		if (m_convMode == CONVERT_NumMarkers)
			num = static_cast<int>(m_convParam);
		else
			num = rng().randGeometric(m_convParam);

		// if conversion reaches end of chromosome, it is an recombination event
		if (num == 0 || num >= ind.lociLeft(index))
			return 0;
		else
			return num;
	} else {
		double len = 0;
		if (m_convMode == CONVERT_TractLength)
			len = m_convParam;
		else
			len = rng().randExponential(len);
		//
		// recombination starts 'before' index so we assume that it happens
		// randomly (uniformly) between this and previous marker
		if (index > 0)
			len -= rng().randUniform01() * ind.lociDist(index - 1, index);
		if (len <= 0. || len >= ind.distLeft(index))
			return 0;
		else
			return ind.lociCovered(index, len);
	}
}


// this function implement how to recombine
// parental chromosomes and set one copy of offspring chromosome
// bt contains the bernulli trailer
void recombinator::recombine(
                             const individual & parent,                                 // one of the parent
                             individual & offspring,                                    // offspring
                             int offPloidy,                                             // which offspring ploidy to fill
                             BernulliTrials & bt,
                             const vectoru & recBeforeLoci,
                             bool setSex)
{
	// use which copy of chromosome
	GenoIterator cp[2], off;

	cp[0] = parent.genoBegin(0);
	cp[1] = parent.genoBegin(1);
	off = offspring.genoBegin(offPloidy);

	// get a new set of values.
	// const BoolResults& bs = bt.trial();
	bt.trial();
	int curCp = bt.trialSucc(recBeforeLoci.size() - 1) ? 0 : 1;
	// the last one does not count, because it determines
	// the initial copy of paternal chromosome
	bt.setTrialSucc(recBeforeLoci.size() - 1, false);

	// algorithm one:
	//
	//  gt: index on chromosomes
	//  gtEnd: total number of loci
	//
	//  at each locus, check if recombine after it, if so
	//  recombine.
	bool withConversion = fcmp_gt(m_convProb, 0);
	if (m_algorithm == 0) {
		// negative means no conversion is pending.
		int convCount = -1;
		size_t gtEnd = recBeforeLoci.back();
		for (size_t gt = 0, bl = 0; gt < gtEnd; ++gt, --convCount) {
			off[gt] = cp[curCp][gt];
			//
			if (convCount == 0) { // conversion ...
				curCp = (curCp + 1) % 2;
				// this is not recorded in m_recCount[bl]
				// no pending conversion
				convCount = -1;
			}
			if (gt + 1 == recBeforeLoci[bl]) {
				// recombination (if convCount == 0, a conversion event is ending)
				if (convCount < 0 && bt.trialSucc(bl)) {
					curCp = (curCp + 1) % 2;
					DBG_DO_(m_recCount[bl]++);
					// if conversion happens
					if (withConversion &&
					    parent.lociLeft(gt) != 1 && // can not be at the end of a chromosome
					    (m_convProb == 1. || rng().randUniform01() < m_convProb)) {
						// convCount will be decreased, until reconversion completes
						// or another recombination happens
						convCount = markersConverted(gt + 1, parent);
						DBG_DO_(m_convSize[convCount]++);
					} else
						// another recombination stops the previous conversion
						convCount = -1;
				}
				++bl;
			}
		}
	} else {
#ifndef BINARYALLELE
		size_t gt = 0, gtEnd = 0;
		size_t pos = bt.probFirstSucc();
		// if there is some recombination
		int convCount = -1;
		size_t convEnd;
		if (pos != BernulliTrials::npos) {
			// first piece
			for (; gt < recBeforeLoci[pos]; ++gt)
				off[gt] = cp[curCp][gt];
			DBG_DO_(m_recCount[pos]++);
			curCp = (curCp + 1) % 2;
			//
			if (withConversion &&
			    parent.lociLeft(gt - 1) != 1 && // can not be at the end of a chromosome
			    (m_convProb == 1. || rng().randUniform01() < m_convProb)) {
				convCount = markersConverted(gt, parent);
				DBG_DO_(m_convSize[convCount]++);
			}
			// next recombination point...
			while ((pos = bt.probNextSucc(pos)) != BernulliTrials::npos) {
				// copy from last to this recombination point, but
				// there might be a conversion event in between
				gtEnd = recBeforeLoci[pos];
				if (convCount > 0) {
					convEnd = gt + convCount;
					if (convEnd < gtEnd) {
						for (; gt < convEnd; ++gt)
							off[gt] = cp[curCp][gt];
						curCp = (curCp + 1) % 2;
					}
					// no pending conversion
					convCount = -1;
				}
				// copy from the end of conversion to this recombination point
				for (; gt < gtEnd; ++gt)
					off[gt] = cp[curCp][gt];
				DBG_DO_(m_recCount[pos]++);
				curCp = (curCp + 1) % 2;
				//
				// conversion event for this recombination event
				if (withConversion &&
				    parent.lociLeft(gt - 1) != 1 && // can not be at the end of a chromosome
				    (m_convProb == 1. || rng().randUniform01() < m_convProb)) {
					// convCount will be decreased, until reconversion completes
					// or another recombination happens
					convCount = markersConverted(gt, parent);
					DBG_DO_(m_convSize[convCount]++);
				}
			}
		}
		gtEnd = recBeforeLoci.back();
		// copy the last piece
		if (convCount > 0) {
			convEnd = gt + convCount;
			if (convEnd < gtEnd) {
				for (; gt < convEnd; ++gt)
					off[gt] = cp[curCp][gt];
				curCp = (curCp + 1) % 2;
			}
		}
		for (; gt < gtEnd; ++gt)
			off[gt] = cp[curCp][gt];
#else
		size_t gt = 0, gtEnd = 0;
		size_t pos = bt.probFirstSucc();
		// if there is some recombination
		int convCount = -1;
		size_t convEnd;
		if (pos != BernulliTrials::npos) {
			// first piece
			gtEnd = recBeforeLoci[pos];
			copyGenotype(cp[curCp] + gt, off + gt, recBeforeLoci[pos] - gt);
			gt = gtEnd;
			DBG_DO_(m_recCount[pos]++);
			curCp = (curCp + 1) % 2;
			if (withConversion &&
			    parent.lociLeft(gt - 1) != 1 && // can not be at the end of a chromosome
			    (m_convProb == 1. || rng().randUniform01() < m_convProb)) {
				convCount = markersConverted(gt, parent);
				DBG_DO_(m_convSize[convCount]++);
			}
			// next recombination point...
			while ((pos = bt.probNextSucc(pos)) != BernulliTrials::npos) {
				gtEnd = recBeforeLoci[pos];
				if (convCount > 0) {
					convEnd = gt + convCount;
					if (convEnd < gtEnd) {
						copyGenotype(cp[curCp] + gt, off + gt, convCount);
						gt = convEnd;
						curCp = (curCp + 1) % 2;
					}
					// no pending conversion
					convCount = -1;
				}
				// copy from the end of conversion to the next recombination point
				copyGenotype(cp[curCp] + gt, off + gt, recBeforeLoci[pos] - gt);
				gt = gtEnd;
				DBG_DO_(m_recCount[pos]++);
				curCp = (curCp + 1) % 2;
				// conversion event for this recombination event
				if (withConversion &&
				    parent.lociLeft(gt - 1) != 1 && // can not be at the end of a chromosome
				    (m_convProb == 1. || rng().randUniform01() < m_convProb)) {
					// convCount will be decreased, until reconversion completes
					// or another recombination happens
					convCount = markersConverted(gt, parent);
					DBG_DO_(m_convSize[convCount]++);
				}
			}
		}
		gtEnd = recBeforeLoci.back();
		// copy the last piece
		if (convCount > 0) {
			convEnd = gt + convCount;
			if (convEnd < gtEnd) {
				copyGenotype(cp[curCp] + gt, off + gt, convCount);
				gt = convEnd;
				curCp = (curCp + 1) % 2;
			}
		}
		copyGenotype(cp[curCp] + gt, off + gt, gtEnd - gt);
#endif
	}
	if (setSex && m_chromX != -1)
		// sex chrom determination
		// if curCp (last chromosome) is X, Female, otherwise Male.
		// Note that for daddy, the last one is arranged XY
		offspring.setSex(curCp == 0 ? Female : Male);
}


// copy the first copy of chromosome from parent to offspring
void recombinator::copyParentalGenotype(const individual & parent,
                                        individual & it,
                                        int ploidy)
{
	GenoIterator par = parent.genoBegin(0);
	GenoIterator off = it.genoBegin(ploidy);

#ifndef BINARYALLELE
	size_t gt = 0;
	size_t gt_end = parent.totNumLoci();
	for (; gt < gt_end; ++gt)
		off[gt] = par[gt];
#else
	copyGenotype(par, off, parent.totNumLoci());
#endif
}


void recombinator::initialize(const population & pop)
{
	m_chromX = pop.chromX();
	m_chromY = pop.chromY();
	m_mitochondrial = pop.mitochondrial();

	// prepare m_bt
	// female
	vectorf vecP;
	// female does not determine sex
	prepareRecRates(pop, m_intensity, m_rate, m_afterLoci,
		m_recBeforeLoci, vecP, Female);

	m_bt.setParameter(vecP, pop.popSize());

	vecP.clear();
	// male case is most complicated.
	double maleIntensity = (m_maleIntensity != -1 || !m_maleRate.empty() )
	                       ? m_maleIntensity : m_intensity;
	vectorf & maleRate = (m_maleIntensity != -1 || !m_maleRate.empty() )
	                     ? m_maleRate : m_rate;
	vectoru & maleAfterLoci = m_maleAfterLoci.empty() ?
	                          m_afterLoci : m_maleAfterLoci;
	// prepare male recombination
	prepareRecRates(pop, maleIntensity, maleRate, maleAfterLoci,
		m_maleRecBeforeLoci, vecP, Male);
	m_maleBt.setParameter(vecP, pop.popSize());
	// choose an algorithm
	// if recombinations are dense. use the first algorithm
	// For example 10 chromoes, regular 0.5*10=5
	// if there are high recombination on chromosomes, ....
	if (std::accumulate(vecP.begin(), vecP.end(), 0.) > pop.numChrom())
		m_algorithm = 0;
	else
		m_algorithm = 1;
	DBG_DO(DBG_RECOMBINATOR, cout << "Algorithm " << m_algorithm << " is being used " << endl);
}


void recombinator::produceOffspring(const individual & parent,
                                    individual & off)
{
	DBG_FAILIF(m_recBeforeLoci.empty(), ValueError,
		"Please use initialize(pop) to set up recombination parameter first.");

	recombine(parent, off, 0, m_bt, m_recBeforeLoci, false);
	recombine(parent, off, 1, m_bt, m_recBeforeLoci, true);
}


void recombinator::produceOffspring(const individual & mom,
                                    const individual & dad, individual & off)
{
	DBG_FAILIF(m_recBeforeLoci.empty(), ValueError,
		"Please use initialize(pop) to set up recombination parameter first.");

	// allows selfing. I.e., if mom or dad is NULL, the other parent will
	// produce both copies of the offspring chromosomes.
	recombine(mom, off, 0, m_bt, m_recBeforeLoci, false);

	if (mom.isHaplodiploid())
		copyParentalGenotype(dad, off, 1);
	else
		// only set sex once for offspring
		recombine(dad, off, 1, m_maleBt, m_maleRecBeforeLoci, true);
}


bool recombinator::applyDuringMating(population & pop,
                                     RawIndIterator offspring,
                                     individual * dad,
                                     individual * mom)
{
	DBG_FAILIF(dad == NULL && mom == NULL, ValueError, "Neither dad or mom is invalid.");

	// call initialize if the signature of pop has been changed.
	baseOperator::applyDuringMating(pop, offspring, dad, mom);

	DBG_FAILIF(m_recBeforeLoci.empty(), ValueError,
		"Uninitialized recombinator");

	if (mom == NULL)
		produceOffspring(*dad, *offspring);
	else if (dad == NULL)
		produceOffspring(*mom, *offspring);
	else
		produceOffspring(*mom, *dad, *offspring);

	return true;
}


}
