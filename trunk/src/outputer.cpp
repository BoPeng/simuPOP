/**
 *  $File: outputer.cpp $
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
#include "outputer.h"

namespace simuPOP {

bool pyOutput::apply(population & pop)
{
	ostream & out = this->getOstream(pop.dict());
	out << m_string;
	this->closeOstream();
	return true;
}


void dumper::displayStructure(const population & pop, ostream & out)
{
	out << "Ploidy: " << pop.ploidy()
	    << " (" << pop.ploidyName() << ")" << endl;
	out << "Chromosomes:\n";
	for (UINT ch = 0; ch < pop.numChrom(); ++ch) {
		out << (ch + 1) << ": " << pop.chromName(ch);
		switch (pop.chromType(ch)) {
		case Autosome:
			out << " (Autosome, ";
			break;
		case ChromosomeX:
			out << " (ChromosomeX, ";
			break;
		case ChromosomeY:
			out << " (ChromosomeY, ";
			break;
		case Customized:
			out << " (Customized, ";
			break;
		default:
			throw ValueError("Wrong chromosome type");
		}
		out << pop.numLoci(ch) << " loci)";
		for (UINT i = 0; i < pop.numLoci(ch); ++i) {
			if (i != 0)
				out << ",";
			if (i % 5 == 0)
				out << "\n ";
			out << " " << pop.locusName(pop.absLocusIndex(ch, i) ) << " ("
			    << pop.locusPos(pop.absLocusIndex(ch, i)) << ")";
		}
		out << endl;
	}
	if (pop.infoSize() != 0) {
		out << "Information fields: " << endl;
		for (UINT i = 0; i < pop.infoSize(); ++i)
			out << pop.infoField(i) << " ";
		out << endl;
	}
	out << "population size: " << pop.popSize();
	out << " (" << pop.numSubPop() << " subpopulations with ";
	for (UINT i = 0, iEnd = pop.numSubPop(); i < iEnd;  ++i) {
		out << (i == 0 ? "" : ", ") << pop.subPopSize(i);
		if (pop.subPopName(i) != UnnamedSubPop)
			out << " (" << pop.subPopName(i) << ")";
	}
	out << " individuals)" << endl;
	out << "Number of ancestral populations: " << pop.ancestralGens() << endl << endl;
}


UINT dumper::displayGenotype(const population & pop, const subPopList & subPops, ostream & out)
{
	UINT count = 0;
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();

	for ( ; sp != spEnd; ++sp) {
		ULONG spSize = pop.subPopSize(*sp);
		out << "Subpopulation " << *sp << " (" << pop.subPopName(*sp) << "), "
		    << toStr(spSize) << " individuals:" << endl;

		if (sp->isVirtual())
			const_cast<population &>(pop).activateVirtualSubPop(*sp, IteratableInds);
		IndIterator ind = const_cast<population &>(pop).indBegin(sp->subPop(), sp->isVirtual() ? IteratableInds : AllInds);
		IndIterator indEnd = const_cast<population &>(pop).indEnd(sp->subPop(), sp->isVirtual() ? IteratableInds : AllInds);

		for ( ; ind != indEnd; ++ind, ++count) {
			out << setw(4) << (& * ind - & * pop.rawIndBegin()) << ": ";
			ind->display(out, m_width, m_loci);
			out << endl;
			if (m_max > 0 && count + 1 >= m_max && count < pop.popSize())
				return count;
		}
	}
	return count;
}


bool dumper::apply(population & pop)
{
	ostream & out = this->getOstream(pop.dict());

	// dump population structure
	if (m_showStructure)
		displayStructure(pop, out);

	subPopList subPops = applicableSubPops();
	if (subPops.empty())
		for (UINT i = 0; i < pop.numSubPop(); ++i)
			subPops.push_back(i);

	if (m_showGenotype) {
		UINT cnt = displayGenotype(pop, subPops, out);
		if (m_max > 0 && cnt == m_max && cnt < pop.popSize())
			out << " ... (" << m_max << " out of " << pop.popSize() << ").\n" << endl;

		int ancGen = m_ancGen;
		// ancGen can be -1
		if (ancGen < 0 || ancGen > static_cast<int>(pop.ancestralGens()))
			ancGen = pop.ancestralGens();
		for (int gen = 1; gen <= ancGen; ++gen) {
			pop.useAncestralGen(gen);
			out << endl << "Ancestry population " << gen << endl;
			UINT cnt = displayGenotype(pop, subPops, out);
			if (m_max > 0 && cnt == m_max && cnt < pop.popSize())
				out << " ... (" << m_max << " out of " << pop.popSize() << ").\n" << endl;

			out << endl;
		}                                                                         // next ancestry
		// IMPORTANT. Reset ancestral pop
		pop.useAncestralGen(0);
	}
	this->closeOstream();
	return true;
}


bool savePopulation::apply(population & pop)
{
	string filename;

	if (m_filename[0] != '!')
		filename = m_filename;
	else {
		Expression filenameParser(m_filename.substr(1));
		filenameParser.setLocalDict(pop.dict());
		filename = filenameParser.valueAsString();
	}
	DBG_DO(DBG_OUTPUTER, cout << "Save to file " << filename << endl);
	pop.save(filename);
	return true;
}


}
