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

namespace simuPOP
{
	bool dumper::apply(population& pop)
	{
		ostream& out = this->getOstream(pop.dict());

		/// dump population structure
		if(! alleleOnly() )
		{
			out << "Ploidy:         \t" << pop.ploidy() << endl;
			out << "Number of chrom:\t" << pop.numChrom() << endl;
			out << "Number of loci: \t";
			for(UINT i=0; i < pop.numChrom(); ++i)
				out << pop.numLoci(i) << " ";
			out << endl;
			out << "Maximum allele state:\t" << pop.maxAllele() << endl;
			out << "Loci positions: " << endl;
			for(UINT ch=0; ch < pop.numChrom(); ++ch)
			{
				cout << "\t" << pop.chromName(ch) << "\t";
				for(UINT i=0; i < pop.numLoci(ch); ++i)
					out << pop.locusPos( pop.absLocusIndex(ch,i) ) << " ";
				out << endl;
			}
			out << "Loci names: " << endl;
			for(UINT ch=0; ch < pop.numChrom(); ++ch)
			{
				out << "\t" << pop.chromName(ch) << "\t";
				for(UINT i=0; i < pop.numLoci(ch); ++i)
					out << pop.locusName( pop.absLocusIndex(ch,i) ) << " ";
				out << endl;
			}
			out << "population size:\t" << pop.popSize() << endl;
			out << "Number of subPop:\t" << pop.numSubPop() << endl;
			out << "Subpop sizes:   \t";
			for(UINT i=0, iEnd = pop.numSubPop(); i < iEnd;  ++i)
				out << pop.subPopSize(i) << " ";
			out << endl;
			out << "Number of ancestral populations:\t" << pop.ancestralDepth() << endl;
		}

		if(! m_infoOnly)
		{
			/// dump all genotypic info
			if(pop.maxAllele() >= 10 && pop.maxAllele() < 100)
				m_width = 2;
			else if( pop.maxAllele() >=100)
				m_width = 3;

			// get individual ranges from subpop
			vectorlu range = m_indRange;
			if(m_indRange.empty())
			{
				if( m_subPop.empty() )			  // all subpop
				{
					for(UINT sp=0; sp < pop.numSubPop();  sp++)
					{
						if( pop.subPopSize(sp) == 0)
							continue;
						range.push_back( pop.subPopBegin(sp));
						range.push_back( pop.subPopEnd(sp));
					}
				}
				else
				{
					for(vectoru::iterator sp=m_subPop.begin();  sp != m_subPop.end();  sp++)
					{
						if( pop.subPopSize(*sp) == 0)
							continue;
						range.push_back( pop.subPopBegin(*sp));
						range.push_back( pop.subPopEnd(*sp));
					}
				}
			}
			out << "individual info: " << endl;
			UINT count = 0;
			for(size_t i=0; i < range.size(); i+=2)
			{
				UINT sp = pop.subPopIndPair(range[i]).first;
				out << "sub population " << sp << ":" << endl;

				for( population::IndIterator ind = pop.indBegin()+range[i];
					ind != pop.indBegin()+range[i+1]; ++ind, ++count)
				{
					out << setw(4) << (ind - pop.indBegin()) << ": ";
					ind->display(out, m_width, m_chrom, m_loci);
					out << endl;

					if(m_max > 0 && count > m_max && count < pop.popSize())
					{
						cout << "population size is " << pop.popSize() << " but dumper() only dump "
							<< m_max << " individuals" << endl
							<< "Use parameter max=-1 to output all individuals." << endl;
						goto done;
					}
				}
			}

			done:
			out << "End of individual info." << endl << endl;

			if(!m_dispAncestry)
			{
				if( pop.ancestralDepth() == 0)
					out << endl << "No ancenstral population recorded." << endl;
				else
					out << endl << "Ignoring " << pop.ancestralDepth() << " ancenstral population(s)." << endl;
			}
			else
			{
				for(size_t i=0; i<pop.ancestralDepth(); ++i)
				{
					pop.useAncestralPop(i+1);
					out << endl << "Ancestry population " << i+1 << endl;

					out << "population size:\t" << pop.popSize() << endl;
					out << "Number of subPop:\t" << pop.numSubPop() << endl;
					out << "Subpop sizes:   \t";
					for(UINT i=0, iEnd = pop.numSubPop(); i < iEnd;  ++i)
						out << pop.subPopSize(i) << " ";
					out << endl;

					out << "individual info: " << endl;

					// get individual ranges from subpop
					vectorlu range = m_indRange;
					if(m_indRange.empty())
					{
												  // all subpop
						if( m_subPop.empty() )
						{
							for(UINT sp=0; sp < pop.numSubPop();  sp++)
							{
								if( pop.subPopSize(sp) == 0)
									continue;
								range.push_back( pop.subPopBegin(sp));
								range.push_back( pop.subPopEnd(sp));
							}
						}
						else
						{
							for(vectoru::iterator sp=m_subPop.begin();  sp != m_subPop.end();  sp++)
							{
								if( pop.subPopSize(*sp) == 0)
									continue;
								range.push_back( pop.subPopBegin(*sp));
								range.push_back( pop.subPopEnd(*sp));
							}
						}
					}
					out << "individual info: " << endl;
					UINT count = 0;
					for(size_t j=0; j < range.size(); j+=2)
					{
						UINT sp = pop.subPopIndPair(range[j]).first;
						out << "sub population " << sp << ":" << endl;

						for( population::IndIterator ind = pop.indBegin()+range[j]; ind != pop.indBegin()+range[j+1]; ++ind)
						{
							out << setw(4) << count++ << ": " ;
							ind->display(out, m_width, m_chrom, m_loci);
							out << endl;

							if(m_max > 0 && count > m_max && count < pop.popSize())
							{
								cout << "population size is " << pop.popSize() << " but dumper() only dump "
									<< m_max << " individuals" << endl
									<< "Use parameter max=-1 to output all individuals." << endl;
								goto doneAnces;
							}
						}
					}

					doneAnces:
					out << "End of individual info." << endl << endl;
				}								  // next ancestry
				// IMPORTANT. Reset ancestral pop
				pop.useAncestralPop(0);
			}									  // dispAncestry
		}
		this->closeOstream();
		return true;
	}

	bool savePopulation::apply(population& pop)
	{
		string filename;
		if( m_filename != "")
			filename = m_filename;
		else
		{
			m_filenameParser.setLocalDict(pop.dict());
			filename = m_filenameParser.valueAsString();
		}
		DBG_DO(DBG_OUTPUTER, cout << "Save to file " << filename << endl);
		pop.savePopulation(filename, m_format, m_compress);
		return true;
	}
}
