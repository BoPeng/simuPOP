/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu                                                        *
*                                                                         *
*   $LastChangedDate$
*   $Rev$
*
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
#include "outputer.h"

namespace simuPOP {
bool dumper::apply(population & pop)
{
	ostream & out = this->getOstream(pop.dict());

	// dump population structure
	if (m_showStructure) {
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
		out << " individuals)";
		out << endl;
		out << "Number of ancestral populations: " << pop.ancestralGens() << endl << endl;
	}

	if (m_showGenotype) {
		// get individual ranges from subpop
		vectorlu range = m_indRange;
		if (m_indRange.empty()) {
			if (m_subPop.empty() ) {                          // all subpop
				for (UINT sp = 0; sp < pop.numSubPop();  sp++) {
					if (pop.subPopSize(sp) == 0)
						continue;
					range.push_back(pop.subPopBegin(sp));
					range.push_back(pop.subPopEnd(sp));
				}
			} else {
				for (vectoru::iterator sp = m_subPop.begin();  sp != m_subPop.end();  sp++) {
					if (pop.subPopSize(*sp) == 0)
						continue;
					range.push_back(pop.subPopBegin(*sp));
					range.push_back(pop.subPopEnd(*sp));
				}
			}
		}
		out << "Genotype of individuals in the present generation:" << endl;
		UINT count = 0;
		for (size_t i = 0; i < range.size(); i += 2) {
			UINT sp = pop.subPopIndPair(range[i]).first;
			out << "Subpopulation " << sp << " (" << pop.subPopName(sp) << "):" << endl;

			for (IndIterator ind = pop.indBegin() + range[i];
			     ind != pop.indBegin() + range[i + 1]; ++ind, ++count) {
				out << setw(4) << (ind - pop.indBegin()) << ": ";
				ind->display(out, m_width, m_chrom, m_loci);
				out << endl;

				if (m_max > 0 && count + 1 >= m_max && count < pop.popSize())
					goto done;
			}
		}

done:
		if (m_max > pop.popSize())
			out << "End of individual genotype.\n" << endl;
		else
			out << "End of individual genotype (" << m_max << " out of "
			    << pop.popSize() << ").\n" << endl;

		int ancGen = m_ancGen;
		// ancGen can be -1
		if (ancGen > static_cast<int>(pop.ancestralGens()))
			ancGen = pop.ancestralGens();
		for (size_t gen = 1; gen <= ancGen; ++gen) {
			pop.useAncestralGen(gen);
			out << endl << "Ancestry population " << gen << endl;
			out << "population size: " << pop.popSize() << endl;
			out << " (" << pop.numSubPop() << " subpopulations with ";
			for (UINT i = 0, iEnd = pop.numSubPop(); i < iEnd;  ++i) {
				out << (i == 0 ? "" : ", ") << pop.subPopSize(i);
				if (pop.subPopName(i) != UnnamedSubPop)
					out << " (" << pop.subPopName(i) << ")";
			}
			out << " individuals)\n";

			// get individual ranges from subpop
			vectorlu range = m_indRange;
			if (m_indRange.empty()) {
				// all subpop
				if (m_subPop.empty() ) {
					for (UINT sp = 0; sp < pop.numSubPop();  sp++) {
						if (pop.subPopSize(sp) == 0)
							continue;
						range.push_back(pop.subPopBegin(sp));
						range.push_back(pop.subPopEnd(sp));
					}
				} else {
					for (vectoru::iterator sp = m_subPop.begin();  sp != m_subPop.end();  sp++) {
						if (pop.subPopSize(*sp) == 0)
							continue;
						range.push_back(pop.subPopBegin(*sp));
						range.push_back(pop.subPopEnd(*sp));
					}
				}
			}
			out << "Genotype info: " << endl;
			UINT count = 0;
			for (size_t j = 0; j < range.size(); j += 2) {
				UINT sp = pop.subPopIndPair(range[j]).first;
				out << "sub population " << sp << " (" << pop.subPopName(sp) << "):" << endl;

				for (IndIterator ind = pop.indBegin() + range[j]; ind != pop.indBegin() + range[j + 1]; ++ind) {
					out << setw(4) << count++ << ": " ;
					ind->display(out, m_width, m_chrom, m_loci);
					out << endl;

					if (m_max > 0 && count > m_max && count < pop.popSize()) {
						cout << "population size is " << pop.popSize() << " but dumper() only dump "
						     << m_max << " individuals" << endl
						     << "Use parameter max=-1 to output all individuals." << endl;
						goto doneAnces;
					}
				}
			}

doneAnces:
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

	if (m_filename != "")
		filename = m_filename;
	else {
		m_filenameParser.setLocalDict(pop.dict());
		filename = m_filenameParser.valueAsString();
	}
	DBG_DO(DBG_OUTPUTER, cout << "Save to file " << filename << endl);
	pop.save(filename);
	return true;
}


}
