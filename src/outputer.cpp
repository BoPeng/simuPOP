/**
 *  $File: outputer.cpp $
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
#include "outputer.h"

namespace simuPOP {

bool PyOutput::apply(Population & pop) const
{
	ostream & out = this->getOstream(pop.dict());

	out << m_string;
	this->closeOstream();
	return true;
}


string PyOutput::describe(bool /* format */) const
{
	return "<simuPOP.PyOutput> write '" + \
	       (m_string.size() > 40 ? m_string.substr(0, 40) + "... " : m_string) + \
	       "' to output";
}


void Dumper::displayStructure(const Population & pop, ostream & out) const
{
	out << "Ploidy: " << pop.ploidy()
	    << " (" << pop.ploidyName() << ")" << endl;
	out << "Chromosomes:\n";
	for (UINT ch = 0; ch < pop.numChrom(); ++ch) {
		out << (ch + 1) << ": " << pop.chromName(ch);
		switch (pop.chromType(ch)) {
		case AUTOSOME:
			out << " (AUTOSOME, ";
			break;
		case CHROMOSOME_X:
			out << " (CHROMOSOME_X, ";
			break;
		case CHROMOSOME_Y:
			out << " (CHROMOSOME_Y, ";
			break;
		case CUSTOMIZED:
			out << " (CUSTOMIZED, ";
			break;
		case MITOCHONDRIAL:
			out << " (MITOCHONDRIAL, ";
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
			out << " " << pop.locusName(pop.absLocusIndex(ch, i)) << " ("
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
	for (size_t i = 0, iEnd = pop.numSubPop(); i < iEnd; ++i) {
		out << (i == 0 ? "" : ", ") << pop.subPopSize(i);
		if (pop.subPopName(i) != UnnamedSubPop)
			out << " (" << pop.subPopName(i) << ")";
	}
	out << " Individuals)" << endl;
	out << "Number of ancestral populations: " << pop.ancestralGens() << endl << endl;
}


UINT Dumper::displayGenotype(const Population & pop, const subPopList & subPops, ostream & out) const
{
	UINT count = 0;
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();

	vectoru infoIdx;
	if (infoSize(&pop) > 0)
		for (size_t i = 0; i < infoSize(); ++i)
			infoIdx.push_back(pop.infoIdx(infoField(i)));

	for ( ; sp != spEnd; ++sp) {
		size_t spSize = pop.subPopSize(*sp);
		out << "SubPopulation " << *sp << " (" << pop.subPopName(*sp) << "), "
		    << spSize << " Individuals:" << endl;

		const_cast<Population &>(pop).activateVirtualSubPop(*sp);
		IndIterator ind = const_cast<Population &>(pop).indIterator(sp->subPop());
		for ( ; ind.valid(); ++ind, ++count) {
			out << setw(4) << (&*ind - &*pop.rawIndBegin()) << ": ";
			ind->display(out, m_width, m_loci, infoIdx);
			out << endl;
			if (m_max > 0 && count + 1 >= m_max && count < pop.popSize())
				break;
		}
		const_cast<Population &>(pop).deactivateVirtualSubPop(sp->subPop());
	}
	return count;
}


bool Dumper::apply(Population & pop) const
{
	ostream & out = this->getOstream(pop.dict());

	// dump population structure
	if (m_showStructure)
		displayStructure(pop, out);

	size_t oldGen = pop.curAncestralGen();
	if (m_showGenotype) {
		vectoru gens = m_ancGens.elems();
		if (m_ancGens.allAvail())
			for (int gen = 0; gen <= pop.ancestralGens(); ++gen)
				gens.push_back(gen);
		else if (m_ancGens.unspecified())
			gens.push_back(pop.curAncestralGen());

		for (unsigned genIdx = 0; genIdx < gens.size(); ++genIdx) {
			pop.useAncestralGen(gens[genIdx]);
			subPopList subPops = applicableSubPops(pop);

			if (gens[genIdx] != 0)
				out << "Ancestral population " << gens[genIdx] << endl;
			UINT cnt = displayGenotype(pop, subPops, out);
			if (m_max > 0 && cnt == m_max && cnt < pop.popSize())
				out << " ... (" << m_max << " out of " << pop.popSize() << ").\n" << endl;

			out << endl;
		}                                                                         // next ancestry
		// IMPORTANT. Reset ancestral pop
		pop.useAncestralGen(oldGen);
	}
	this->closeOstream();
	return true;
}


string SavePopulation::describe(bool /* format */) const
{
	return "<simuPOP.SavePopulation> save population to file " + m_filename;
}


bool SavePopulation::apply(Population & pop) const
{
	if (m_filename.empty())
		return true;

	string filename;

	if (m_filename[0] != '!')
		filename = m_filename;
	else {
		Expression filenameParser(m_filename.substr(1));
		filenameParser.setLocalDict(pop.dict());
		filename = filenameParser.valueAsString();
	}
	DBG_DO(DBG_POPULATION, cerr << "Save to file " << filename << endl);
	pop.save(filename);
	return true;
}


}
