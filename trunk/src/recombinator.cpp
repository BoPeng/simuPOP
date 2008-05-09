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

namespace simuPOP {
void recombinator::prepareRecRates(population & pop,
                                   double intensity,
                                   vectorf rate,
                                   vectoru afterLoci,                   //
                                   bool sexChrom,                       // whether or not recombine the last chromosome
                                   vectoru & recBeforeLoci,             // return before loci vector
                                   vectorf & vecP)                      // return recombination rate
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

	bool useLociDist;
	if (rate.empty() )  // actually use intensity
		useLociDist = true;
	else
		useLociDist = false;

	// first, we create a recBeforeLoci vector
	// and a recombination rate vector
	if (afterLoci.empty()) {
		size_t vecSize = pop.totNumLoci();
		if (sexChrom)
			vecSize -= pop.numLoci(pop.numChrom() - 1) - 1;

		recBeforeLoci.resize(vecSize);
		for (size_t i = 0; i < vecSize; ++i)
			recBeforeLoci[i] = i + 1;
		// if sex chrom, the last one is the end of geno
		// instead of the beginning of the last chromosome
		recBeforeLoci[vecSize - 1] = pop.totNumLoci();

		vecP.resize(vecSize);

		int index = 0;
		UINT chEnd = pop.numChrom();

		for (UINT ch = 0; ch < chEnd; ++ch) {
			// do we recombine the last chromosome?
			if (ch == chEnd - 1 && sexChrom) {
				// the last one is used to determine the choice of
				// the first chromosome copy
				vecP[index] = 0.5;
				break;
			}

			/// get loci distance * rate and then recombinant points
			for (UINT loc = 0, locEnd = pop.numLoci(ch) - 1;
			     loc < locEnd; ++loc) {
				if (useLociDist)
					vecP[index] = (pop.locusPos(index + 1) - pop.locusPos(index)) * intensity;
				else
					vecP[index] = rate[0];

				DBG_WARNING(fcmp_gt(vecP[index], 0.5),
					"Recombination rate after marker " + toStr(index) + " is out of range ("
					+ toStr(vecP[index]) + " ) so it is set to 0.5. This may happen \n"
					"when you use recombination intensity instead of rate, and your loci \n"
					"distance is too high.)");
				if (fcmp_gt(vecP[index], 0.5))
					vecP[index] = 0.5;
				index++;
			}
			// 'recombine' after each chromosome.
			vecP[index++] = .5;
		}

		// FIXME: remove zero indexes for efficiency purpose, but there is no
		// real need for this if atLoci is not specified.

		DBG_DO(DBG_RECOMBINATOR, cout << "Use all Loci. With rates "
			                          << vecP << " before " << recBeforeLoci << endl);
	} else {                                                                          // afterLoci not empty
		DBG_FAILIF(rate.size() > 1 && rate.size() != afterLoci.size(), SystemError,
			"If an array is given, rates and afterLoci should have the same length");

		vectoru::iterator pos;
		vecP.clear();
		recBeforeLoci.clear();

		UINT index = 0;

		for (UINT ch = 0, chEnd = pop.numChrom(); ch < chEnd; ++ch) {
			if (ch == chEnd - 1 && sexChrom) {
				vecP.push_back(0.5);
				recBeforeLoci.push_back(pop.totNumLoci());
				break;
			}

			/// get loci distance * rate and then recombinant points
			for (UINT loc = 0, locEnd = pop.numLoci(ch) - 1;
			     loc < locEnd; ++loc) {
				// if this locus will be recombined.
				pos = find(afterLoci.begin(), afterLoci.end(), index);
				if (pos != afterLoci.end()) {
					if (useLociDist) {
						if (intensity > 0) {              // igore zero rate
							vecP.push_back( (pop.locusPos(index + 1) - pop.locusPos(index)) * intensity);
							recBeforeLoci.push_back(index + 1);
						}
					} else if (rate.size() == 1 && !useLociDist) {
						if (rate[0] > 0) {                // ignore zero rate
							vecP.push_back(rate[0]);
							recBeforeLoci.push_back(index + 1);
						}
					} else {
						// ignore zero rate
						if (rate[pos - afterLoci.begin()] > 0) {
							vecP.push_back(rate[pos - afterLoci.begin()]);
							recBeforeLoci.push_back(index + 1);
						}
					}

					DBG_ASSERT(fcmp_ge(vecP[vecP.size() - 1], 0) && fcmp_le(vecP[vecP.size() - 1], 1),
						ValueError,
						"Recombination rate should be in [0,1]. (Maybe your loci distance is too high.)");
				}
				index++;
			}
			vecP.push_back(.5);

			DBG_WARNING(find(afterLoci.begin(), afterLoci.end(), index) != afterLoci.end(), 
				"Specified recombination rate for the last locus on a chromosome is discarded.");

			// add between chromosomes....
			index++;
			recBeforeLoci.push_back(index);
		}

		// for debug purpose, check the original list:
		for (size_t i = 0; i < afterLoci.size(); ++i) {
			DBG_FAILIF(find(recBeforeLoci.begin(), recBeforeLoci.end(), afterLoci[i] + 1) == recBeforeLoci.end(),
				ValueError, "Specified locus " + toStr(afterLoci[i]) +
				" is not used in recombinator. Is it valid?");
		}

		DBG_ASSERT(vecP.size() == recBeforeLoci.size(), SystemError,
			"Rate and before loci should have the same length.");

		DBG_ASSERT(recBeforeLoci.back() == pop.totNumLoci(),
			SystemError, "The last beforeLoci elem should be total number of loci.");

		DBG_ASSERT(vecP.back() == .5, SystemError,
			"The last elem of rate should be half.");

		DBG_DO(DBG_RECOMBINATOR, cout << "Specify after Loci. With rates "
			                          << vecP << " before " << recBeforeLoci << endl);
	}

	/// initialize recombination counter,
	/// This will count recombination events after
	/// each locus.
	DBG_DO_(m_recCount.resize(pop.totNumLoci(), 0));
	return;
}


int recombinator::markersConverted(size_t index, individual * ind)
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
		if (num == 0 || num >= ind->lociLeft(index))
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
			len -= rng().randUniform01() * ind->lociDist(index - 1, index);
		if (len <= 0. || len >= ind->distLeft(index))
			return 0;
		else
			return ind->lociCovered(index, len);
	}
}


// this function implement how to recombine
// parental chromosomes and set one copy of offspring chromosome
// bt contains the bernulli trailer
void recombinator::recombine(
                             individual * parent,                               // one of the parent
                             RawIndIterator & offspring,                        // offspring
                             int offPloidy,                                     // which offspring ploidy to fill
                             BernulliTrials & bt,
                             const vectoru & recBeforeLoci,
                             bool setSex)
{
	// use which copy of chromosome
	GenoIterator cp[2], off;

	cp[0] = parent->genoBegin(0);
	cp[1] = parent->genoBegin(1);
	off = offspring->genoBegin(offPloidy);

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
					    parent->lociLeft(gt) != 1 && // can not be at the end of a chromosome
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
			    parent->lociLeft(gt - 1) != 1 && // can not be at the end of a chromosome
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
				    parent->lociLeft(gt - 1) != 1 && // can not be at the end of a chromosome
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
			    parent->lociLeft(gt - 1) != 1 && // can not be at the end of a chromosome
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
				    parent->lociLeft(gt - 1) != 1 && // can not be at the end of a chromosome
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
	if (setSex && m_hasSexChrom)
		// sex chrom determination
		// if curCp (last chromosome) is X, Female, otherwise Male.
		// Note that for daddy, the last one is arranged XY
		offspring->setSex(curCp == 0 ? Female : Male);
}


// copy the first copy of chromosome from parent to offspring
void recombinator::copyParentalGenotype(individual * parent,
                                        RawIndIterator & it,
                                        int ploidy)
{
	GenoIterator par = parent->genoBegin(0);
	GenoIterator off = it->genoBegin(ploidy);

#ifndef BINARYALLELE
	size_t gt = 0;
	size_t gt_end = parent->totNumLoci();
	for (; gt < gt_end; ++gt)
		off[gt] = par[gt];
#else
	copyGenotype(par, off, parent->totNumLoci());
#endif
}


bool recombinator::applyDuringMating(population & pop,
                                     RawIndIterator offspring,
                                     individual * dad,
                                     individual * mom)
{
	DBG_FAILIF(dad == NULL && mom == NULL, ValueError, "Neither dad or mom is invalid.");

	// first time setup
	if (m_recBeforeLoci.empty()) {
		// prepare m_bt
		// female
		vectorf vecP;
		// female does not determine sex
		prepareRecRates(pop, m_intensity, m_rate, m_afterLoci,
			false, m_recBeforeLoci, vecP);

		m_bt.setParameter(vecP, pop.popSize());

		vecP.clear();
		// male case is most complicated.
		m_hasSexChrom = pop.sexChrom() ? true : false;
		double maleIntensity = (m_maleIntensity != -1 || !m_maleRate.empty() )
		                       ? m_maleIntensity : m_intensity;
		vectorf & maleRate = (m_maleIntensity != -1 || !m_maleRate.empty() )
		                     ? m_maleRate : m_rate;
		vectoru & maleAfterLoci = m_maleAfterLoci.empty() ?
		                          m_afterLoci : m_maleAfterLoci;
		// prepare male recombination
		prepareRecRates(pop, maleIntensity, maleRate, maleAfterLoci,
			m_hasSexChrom, m_maleRecBeforeLoci, vecP);
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

	// allows selfing. I.e., if mom or dad is NULL, the other parent will
	// produce both copies of the offspring chromosomes.
	if (mom != NULL)
		recombine(mom, offspring, 0, m_bt, m_recBeforeLoci, false);
	else
		recombine(dad, offspring, 0, m_bt, m_recBeforeLoci, false);

	if (pop.haplodiploid()) {
		DBG_FAILIF(dad == NULL, ValueError,
			"Invalid male in haplodiploid mode");
		copyParentalGenotype(dad, offspring, 1);
		return true;
	}

	if (dad != NULL)
		// only set sex once for offspring
		recombine(dad, offspring, 1, m_maleBt, m_maleRecBeforeLoci, true);
	else
		recombine(mom, offspring, 1, m_maleBt, m_maleRecBeforeLoci, true);
	return true;
}


}
