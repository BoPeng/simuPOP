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

void cloneGenoTransmitter::initialize(const population & pop)
{
	m_hasCustomizedChroms = !pop.customizedChroms().empty();
	if (m_hasCustomizedChroms) {
		for (UINT ch = 0; ch < pop.numChrom(); ++ch)
			if (pop.chromType(ch) == Customized)
				m_lociToCopy.push_back(0);
			else
				m_lociToCopy.push_back(pop.numLoci(ch));
	}
	m_ploidy = pop.ploidy();
}


bool cloneGenoTransmitter::applyDuringMating(population & pop,
                                             RawIndIterator offspring,
                                             individual * dad,
                                             individual * mom)
{
	DBG_FAILIF(dad == NULL && mom == NULL, ValueError,
		"Both parents are invalid");

	// call initialize if needed.
	baseOperator::applyDuringMating(pop, offspring, dad, mom);

	individual * parent = mom != NULL ? mom : dad;

	// troublesome ...
	if (m_hasCustomizedChroms) {
		for (UINT p = 0; p != m_ploidy; ++p) {
			for (UINT ch = 0; ch < pop.numChrom(); ++ch) {
				if (m_lociToCopy[ch] == 0)
					continue;
				GenoIterator par = parent->genoBegin(p, ch);
				GenoIterator off = offspring->genoBegin(p, ch);
#ifdef BINARYALLELE
				copyGenotype(par, off, m_lociToCopy[ch]);
#else
				GenoIterator par_end = parent->genoEnd(p, ch);
				copy(par, par_end, off);
#endif
			}
		}
	} else {             // easy
#ifdef BINARYALLELE
		copyGenotype(parent->genoBegin(), offspring->genoBegin(),
			offspring->genoSize());
#else
		copy(parent->genoBegin(), parent->genoEnd(), offspring->genoBegin());
#endif
	}
	return true;
}


void mendelianGenoTransmitter::initialize(const population & pop)
{
	vectorf prob(2 * pop.numChrom(), 0.5);

	m_bt.setParameter(prob, pop.popSize());
	m_chIdx = pop.chromIndex();

	m_hasCustomizedChroms = !pop.customizedChroms().empty();
	for (UINT ch = 0; ch < pop.numChrom(); ++ch) {
		if (pop.chromType(ch) == Customized)
			m_lociToCopy.push_back(0);
		else
			m_lociToCopy.push_back(pop.numLoci(ch));
	}
	m_chromX = pop.chromX();
	m_chromY = pop.chromY();
	m_numChrom = pop.numChrom();
}


void mendelianGenoTransmitter::formOffspringGenotype(individual * parent,
                                                     RawIndIterator & it, int ploidy)
{
	// current parental ploidy (copy from which chromosome copy)
	int parPloidy = 0;
	// pointer to parental, and offspring chromosome copies
	GenoIterator par[2];
	GenoIterator off;

	//
	par[0] = parent->genoBegin(0);
	par[1] = parent->genoBegin(1);
	off = it->genoBegin(ploidy);
	//
	int btShift = ploidy * m_numChrom;
#ifndef BINARYALLELE
	// the easy way to copy things.
	for (UINT ch = 0; ch < m_numChrom; ++ch) {
		// customized chromosome?
		if (m_lociToCopy[ch] == 0 ||
		    (ploidy == 0 && ch == m_chromY))    // maternal, Y chromosome
			continue;
		if (ploidy == 1 && ch == m_chromX) {
			if (it->sex() == Male)
				continue;
			else
				parPloidy = 0;
		} else if (ploidy == 1 && ch == m_chromY) {
			if (it->sex() == Male)
				parPloidy = 1;          // copy chrom Y from second ploidy
			else
				continue;
		} else
			parPloidy = m_bt.trialSucc(ch + btShift);
		//
		for (size_t gt = m_chIdx[ch]; gt < m_chIdx[ch + 1]; ++gt)
			off[gt] = par[parPloidy][gt];
	}
#else
	// for the simple case, use faster algorithm
	if (m_chromX < 0 && m_chromY < 0 && !m_hasCustomizedChroms()) {
		//
		// 1. try to copy in blocks,
		// 2. if two chromosomes can be copied together, copy together
		// 3. if length is short, using the old method.
		//
		size_t parBegin = 0;
		size_t parEnd = 0;
		// first chromosome
		parPloidy = m_bt.trialSucc(btShift);
		//
		int nextParPloidy = 0;
		bool copyPar;
		for (UINT ch = 0; ch < chEnd; ++ch) {
			// if it is the last chromosome, copy anyway
			if (ch == chEnd - 1)
				copyPar = true;
			else {                                                                 // is there a different chromosome?
				nextParPloidy = m_bt.trialSucc(ch + 1 + btShift);
				copyPar = parPloidy != nextParPloidy;
			}
			if (copyPar) {
				// end of this chromosome, is the beginning of the next
				parEnd = m_chIdx[ch + 1];
				size_t length = parEnd - parBegin;
				//
				// the easiest case, try to get some speed up...
				if (length == 1)
					off[parBegin] = par[parPloidy][parBegin];
				else
					copyGenotype(par[parPloidy] + parBegin, off + parBegin, length);
				//
				if (ch != chEnd - 1)
					parPloidy = nextParPloidy;
				parBegin = parEnd;
			}
		}
	} else {    // use the less efficient algorithm
		for (UINT ch = 0; ch < m_numChrom; ++ch) {
			// customized chromosome?
			if (m_lociToCopy[ch] == 0 ||
			    (ploidy == 0 && ch == m_chromY))    // maternal, Y chromosome
				continue;
			if (ploidy == 1 && ch == m_chromX) {
				if (it->sex() == Male)
					continue;
				else
					parPloidy = 0;
			} else if (ploidy == 1 && ch == m_chromY) {
				if (it->sex() == Male)
					parPloidy = 1;          // copy chrom Y from second ploidy
				else
					continue;
			} else
				parPloidy = m_bt.trialSucc(ch + btShift);
			//
			copyGenotype(par[parPloidy] + gt, off + gt, m_lociToCopy[ch]);
		}
	}
#endif
}


bool mendelianGenoTransmitter::applyDuringMating(population & pop,
                                                 RawIndIterator offspring, individual * dad, individual * mom)
{
	DBG_FAILIF(mom == NULL || dad == NULL, ValueError,
		"Mendelian offspring generator requires two valid parents");

	// call initialize if needed.
	baseOperator::applyDuringMating(pop, offspring, dad, mom);

	// m_bt 's width is 2*numChrom() and can be used for
	// the next two functions.
	m_bt.trial();
	formOffspringGenotype(mom, offspring, 0);
	formOffspringGenotype(dad, offspring, 1);
	return true;
}


bool selfingGenoTransmitter::applyDuringMating(population & pop,
                                               RawIndIterator offspring, individual * dad, individual * mom)
{
	//
	DBG_FAILIF(mom == NULL || dad == NULL, ValueError,
		"Mendelian offspring generator requires two valid parents");

	// call initialize if needed.
	baseOperator::applyDuringMating(pop, offspring, dad, mom);

	individual * parent = mom != NULL ? mom : dad;

	// m_bt 's width is 2*numChrom() and can be used for
	// the next two functions.
	m_bt.trial();
	// use the same parent to produce two copies of chromosomes
	formOffspringGenotype(parent, offspring, 0);
	formOffspringGenotype(parent, offspring, 1);
	return true;
}


bool haplodiploidGenoTransmitter::applyDuringMating(population & pop,
                                                    RawIndIterator offspring, individual * dad, individual * mom)
{
	DBG_FAILIF(dad == NULL || mom == NULL, ValueError,
		"haplodiploid offspring generator: one of the parents is invalid.");

	// call initialize if needed.
	baseOperator::applyDuringMating(pop, offspring, dad, mom);

	// m_bt 's width is 2*numChrom() and can be used for
	// the next two functions.
	m_bt.trial();
	// mom generate the first...
	formOffspringGenotype(mom, offspring, 0);

	//
	if (offspring->sex() == Female) {
		// paternal
		GenoIterator par = dad->genoBegin(0);
		GenoIterator off = offspring->genoBegin(1);

		//
#ifndef BINARYALLELE
		size_t gt = 0;
		size_t gt_end = dad->totNumLoci();
		for (; gt < gt_end; ++gt)
			off[gt] = par[gt];
#else
		copyGenotype(par, off, dad->totNumLoci());
#endif
	}
	return true;
}


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

		if (std::find(m_customized.begin(), m_customized.end(), ch) != m_customized.end() ||             // customized does not recombine
		    (sex == Male && (static_cast<int>(ch) == m_chromX || static_cast<int>(ch) == m_chromY)) ||
		    (sex == Female && static_cast<int>(ch) == m_chromY)) {
			DBG_DO(DBG_RECOMBINATOR, cout << "Ignoring chromosome " << ch
				                          << " Customized: " << m_customized << " Chrom X: " << m_chromX
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
			if (convCount == 0) {             // conversion ...
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
					    parent.lociLeft(gt) != 1 &&             // can not be at the end of a chromosome
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
			    parent.lociLeft(gt - 1) != 1 &&             // can not be at the end of a chromosome
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
				    parent.lociLeft(gt - 1) != 1 &&             // can not be at the end of a chromosome
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
			    parent.lociLeft(gt - 1) != 1 &&             // can not be at the end of a chromosome
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
				    parent.lociLeft(gt - 1) != 1 &&             // can not be at the end of a chromosome
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
	m_customized = pop.customizedChroms();

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
