/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
 *                                                                         *
 *   $LastChangedDate: 2006-02-21 15:27:25 -0600 (Tue, 21 Feb 2006) $
 *   $Rev: 191 $
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

#include "initializer.h"

namespace simuPOP
{

	void initializer::setRanges(population& pop)
	{
		m_ranges = m_indRange;

		if( m_ranges.empty())
		{
			if( m_subPop.empty())
			{
				m_ranges.resize( pop.numSubPop() );
				for(size_t i=0; i < pop.numSubPop(); ++i)
				{
					m_ranges[i].resize(2, pop.subPopBegin(i));
					m_ranges[i][1] = pop.subPopEnd(i);
				}
			}
			else
			{
				m_ranges.resize( m_subPop.size());
				for(size_t i = 0; i < m_subPop.size(); ++i)
				{
					m_ranges[i].resize(2, pop.subPopBegin(m_subPop[i]));
					m_ranges[i][1] = pop.subPopEnd(m_subPop[i]);
				}
			}
		}
	}

	Sex initializer::nextSex()
	{
		DBG_ASSERT(*m_sexItr == int(Male) || *m_sexItr == int(Female),
			ValueError, "sex must be array of Male or Female. "
			+ toStr(*m_sexItr) + " given.");
		Sex s = *m_sexItr++ == 1?Male:Female;
		if(m_sexItr == m_sex.end())
			m_sexItr = m_sex.begin();
		return s;
	}

	bool initByFreq::apply(population& pop)
	{
		/// initialize m_ranges
		setRanges(pop);
#ifdef SIMUMPI
		// still run without marker, go trigger errors etc.
		bool noMarker = false;
		// for mpi, there is no data in the first node
		if (mpiRank() == 0)
			noMarker = true;
#endif

		initSexIter();

		DBG_FAILIF( m_alleleFreq.size() > 1 && m_alleleFreq.size() != m_ranges.size(),
			ValueError, "Ranges and values should have the same length");

#ifdef SIMUMPI
		// if atLoci is given, fit to current node
		if (!m_atLoci.empty())
		{
			vectori localLoci;
			for(size_t loc=0; loc < m_atLoci.size(); ++loc)
				if (pop.rankOfLocus(m_atLoci[loc]) == mpiRank())
					localLoci.push_back(m_atLoci[loc]);
			if(localLoci.empty())
				noMarker = true;
		}
#endif
		for(size_t rg = 0; rg < m_ranges.size(); ++rg)
		{
			vectorf& alleleFreq = m_alleleFreq.size() == 1 ? m_alleleFreq[0] : m_alleleFreq[rg];

			ULONG left = m_ranges[rg][0], right = m_ranges[rg][1];

			DBG_FAILIF( left > pop.popSize() || right > pop.popSize() || left > right ,
				ValueError, "Invaid range boundary: " + toStr(left) + " - " + toStr(right-1));

			// Weightedsampler ws(rng(), incFreq);
			Weightedsampler ws(rng(), alleleFreq);

			DBG_ASSERT( fcmp_eq(std::accumulate(alleleFreq.begin(), alleleFreq.end(), 0.), 1),
				SystemError, "Allele frequecies should add up to one.");

			pop.adjustGenoPosition(true);
			if(m_identicalInds)
			{
				if(m_atLoci.empty())			  // to all loci
				{
					if( m_atPloidy==-1)			  // all chromosomes
					{
						// for MPI mode, indGenoBegin and indGenoEnd
						// refers to local pieces.
						ws.get(pop.indGenoBegin(left), pop.indGenoEnd(left));

						for(ULONG ind=left+1; ind != right; ++ind)
							copy(pop.indGenoBegin(left), pop.indGenoEnd(left), pop.indGenoBegin(ind));
					}
					else						  // only initialize one set of chromosome
					{
						// for MPI mode, indGenoBegin and indGenoEnd
						// refers to local pieces.
						ws.get(pop.ind(left).genoBegin(m_atPloidy),
							pop.ind(left).genoEnd(m_atPloidy));

						for(ULONG ind=left+1; ind != right; ++ind)
							copy(pop.ind(left).genoBegin(m_atPloidy),
								pop.ind(left).genoEnd(m_atPloidy),
								pop.ind(ind).genoBegin(m_atPloidy));
					}
				}
				else
				{
					// initislize locus by locus
#ifdef SIMUMPI					
					if(!noMarker)
					{
#endif						
						for(vectoru::iterator locus=m_atLoci.begin(); locus != m_atLoci.end(); ++locus)
						{
#ifdef SIMUMPI							
							if(pop.rankOfLocus(*locus) != mpiRank())
								continue;
#endif								
							// all chromosomes
							if(m_atPloidy == -1)
							{
								ws.get(pop.alleleBegin(*locus, false) + left * pop.ploidy(),
									pop.alleleBegin(*locus, false) + (left+1) * pop.ploidy() );

								for(ULONG ind=left+1; ind < right; ++ind)
									copy(pop.alleleBegin(*locus, false) + left*pop.ploidy(),
										pop.alleleBegin(*locus, false) + (left+1)*pop.ploidy(),
										pop.alleleBegin(*locus, false) + ind*pop.ploidy());
							}
							else				  // only one of the copies (do it one by one?)
							{
								UINT a = ws.get();
								for(ULONG ind=left; ind != right; ++ind)
									pop.ind(ind).setAllele(a, *locus, m_atPloidy);
							}
						}
#ifdef SIMUMPI						
					}
#endif					
				}
			}
			else								  // not idential individuals
			{
				if( m_atLoci.empty())			  // at all loci
				{
					if( m_atPloidy == -1)
					{
#ifdef SIMUMPI
						ws.get(pop.genoBegin(false)+left*pop.localGenoSize(),
							pop.genoBegin(false)+right*pop.localGenoSize());
#else
						ws.get(pop.genoBegin(false)+left*pop.genoSize(),
							pop.genoBegin(false)+right*pop.genoSize());
#endif
					}
					else						  // for only one ploidy
					{
						for(ULONG ind=left; ind != right; ++ind)
							ws.get( pop.ind(ind).genoBegin(m_atPloidy),
								pop.ind(ind).genoEnd(m_atPloidy));
					}
				}
				else							  // at certain loci
				{
#ifdef SIMUMPI
					if(!noMarker)
					{
#endif						
						if( m_atPloidy == -1)
						{
							for(vectoru::iterator locus=m_atLoci.begin(); locus != m_atLoci.end(); ++locus)
								ws.get(pop.alleleBegin(*locus, false) + left * pop.ploidy(),
									pop.alleleBegin(*locus, false) + right * pop.ploidy() );
						}
						else					  // for only one ploidy
						{
							for(vectoru::iterator locus=m_atLoci.begin(); locus != m_atLoci.end(); ++locus)
							{
#ifdef SIMUMPI
								if(pop.rankOfLocus(*locus) != mpiRank())
									continue;
#endif
								for(ULONG ind=left; ind != right; ++ind)
									pop.ind(ind).setAllele(ws.get(), *locus, m_atPloidy);
							}
						}
#ifdef SIMUMPI						
					}
#endif					
				}
			}
			// initialize sex
			// call randUnif once for each individual
			// (initialize allele need to call randUnif for each locus
			if(m_sex.empty())
			{
				for ( ULONG ind = left; ind != right; ++ind)
				{
					if( rng().randUniform01() < m_maleFreq )
						pop.ind(ind).setSex( Male );
					else
						pop.ind(ind).setSex( Female );
				}
			}
			else								  // use provided sex array
			{
				for ( ULONG ind = left; ind != right; ++ind)
				{
					pop.ind(ind).setSex( nextSex());
				}
			}									  // set sex
		}										  // range
		return true;
	}

	bool initByValue::apply(population& pop)
	{
#ifdef SIMUMPI
		// still run without marker, go trigger errors etc.
		bool noMarker = false;
		// for mpi, there is no data in the first node
		if (mpiRank() == 0)
			noMarker = true;
#endif

		initSexIter();

		for(size_t i = 0; i < m_value.size(); ++i)
		{
			// extend m_value
			if(m_proportion.empty())
			{
				if(m_atLoci.empty() && m_value[i].size() == pop.totNumLoci())
				{
#ifdef SIMUMPI
					if( m_atPloidy==-1)
					{
						// shift pieces to the beginning
						vectori tmpValue(pop.localNumLoci()*pop.ploidy());
						for(UINT p=0; p < pop.ploidy(); ++p)
							copy(m_value[i].begin() + pop.beginLocus(),
								m_value[i].begin() + pop.endLocus(),
								tmpValue.begin() + pop.localNumLoci()*p);
						m_value[i].swap(tmpValue);
					}
					else
					{
						vectori tmpValue(pop.localNumLoci());
						copy(m_value[i].begin() + pop.beginLocus(),
							m_value[i].begin() + pop.endLocus(),
							tmpValue.begin());
						m_value[i].swap(tmpValue);
					}
#else
					if( m_atPloidy==-1)
					{
						m_value[i].resize( pop.totNumLoci() * pop.ploidy());
						for(UINT p=1; p < pop.ploidy(); ++p)
							copy( m_value[i].begin(),
								m_value[i].begin() + pop.totNumLoci(),
								m_value[i].begin() + pop.totNumLoci()*p);
					}
#endif
				}
#ifdef SIMUMPI
				if(m_atLoci.empty() && m_value[i].size() == pop.totNumLoci()*pop.ploidy())
				{
					if(m_atPloidy==-1)
					{
						vectori tmpValue(pop.localNumLoci()*pop.ploidy());
						for(UINT p=0; p < pop.ploidy(); ++p)
							copy(m_value[i].begin() + pop.beginLocus() + pop.totNumLoci()*p,
								m_value[i].begin() + pop.endLocus() + pop.totNumLoci()*p,
								tmpValue.begin() + pop.localNumLoci()*p);
						m_value[i].swap(tmpValue);
					}
				}
#endif
				if( !m_atLoci.empty() &&
					(m_value[i].size() == m_atLoci.size() ||
					m_value[i].size() == m_atLoci.size() * pop.ploidy())
					)
				{
#ifdef SIMUMPI
					// pick at loci
					vectoru tmpLoci;
					vectori tmpValue;
					for(size_t loc=0; loc < m_atLoci.size(); ++loc)
					{
						if (pop.rankOfLocus(m_atLoci[loc]) == mpiRank())
						{
							tmpLoci.push_back(m_atLoci[loc] - pop.beginLocus());
							tmpValue.push_back(m_value[i][loc]);
						}
					}
					m_atLoci.swap(tmpLoci);
					if (m_atLoci.empty())
						noMarker = true;
					if( m_atPloidy==-1)
					{
						tmpValue.resize(tmpLoci.size() * pop.ploidy());
						for(UINT p=1; p < pop.ploidy(); ++p)
							copy(tmpValue.begin(),
								tmpValue.begin() + tmpLoci.size(),
								tmpValue.begin() + tmpLoci.size()*p);
						m_value[i].swap(tmpValue);
					}
#else
					if( m_atPloidy==-1)
					{
						m_value[i].resize( m_atLoci.size() * pop.ploidy());
						for(UINT p=1; p < pop.ploidy(); ++p)
							copy( m_value[i].begin(),
								m_value[i].begin() + m_atLoci.size(),
								m_value[i].begin() + m_atLoci.size()*p);
					}
#endif
				}								  // case 2
			}									  // no proportion
		}										  // check m_value[i]
#ifdef SIMUMPI
		if(!m_proportion.empty() && !m_atLoci.empty())
		{
			for(size_t i = 0; i < m_value.size(); ++i)
			{
				vectori tmpValue;
				for(size_t loc=0; loc < m_atLoci.size(); ++loc)
					if (pop.rankOfLocus(m_atLoci[loc]) == mpiRank())
						tmpValue.push_back(m_value[i][loc]);
				m_value[i].swap(tmpValue);
			}
			vectoru tmpLoci;
			for(size_t loc=0; loc < m_atLoci.size(); ++loc)
				if (pop.rankOfLocus(m_atLoci[loc]) == mpiRank())
					tmpLoci.push_back(m_atLoci[loc] - pop.beginLocus());
			m_atLoci.swap(tmpLoci);
			if (m_atLoci.empty())
				noMarker = true;
		}
#endif

		DBG_DO(DBG_INITIALIZER, cout << "Size of src is " << m_value[0].size() << endl);

#ifndef OPTIMIZED
		UINT gSz= m_value[0].size();
		for(size_t v=1; v< m_value.size(); ++v)
			if( m_value[v].size() != gSz)
				throw ValueError("Given values should have the same length (either one copy of chromosomes or the whole genotype.");
#endif

		setRanges(pop);

		DBG_FAILIF( m_proportion.empty() && m_value.size() > 1
			&& m_value.size() != m_ranges.size(),
			ValueError, "Ranges and values should have the same length");

		if( m_proportion.empty() )
		{
			// for each range
			for(size_t rg = 0; rg < m_ranges.size(); ++rg)
			{
				// we have left and right
				ULONG left = m_ranges[rg][0], right = m_ranges[rg][1];

				DBG_FAILIF( left > pop.popSize() || right > pop.popSize() || left > right ,
					ValueError, "Invaid m_ranges boundary: " + toStr(left) + " - " + toStr(right));

				vectori& src = m_value.size() > 1 ? m_value[rg] : m_value[0];
				size_t srcSz = src.size();
				size_t lociSz = m_atLoci.size();
#ifdef SIMUMPI
				size_t totNumLoci = pop.localNumLoci();
#else
				size_t totNumLoci = pop.totNumLoci();
#endif

				for (ULONG ind = left; ind < right; ++ind)
				{
#ifdef SIMUMPI
					if(m_atLoci.empty() && !noMarker)
#else
						if(m_atLoci.empty())
#endif
					{
						if(m_atPloidy==-1)	  // all copies of chromosome
						{
#ifdef SIMUMPI
							DBG_ASSERT(src.size() == pop.localGenoSize(), ValueError,
								"Rank " + toStr(mpiRank()) + " Length of value does not match geno size");
#else
							DBG_ASSERT(src.size() == pop.genoSize(), ValueError,
								"Length of value does not match geno size");
#endif
							copy(src.begin(), src.end(), pop.indGenoBegin(ind));
						}
						else					  // one of the copied.
						{						  /// fixme: check length of src?
#ifdef SIMUMPI
							DBG_ASSERT(src.size() == pop.localNumLoci(), ValueError,
								"Rank " + toStr(mpiRank()) + "Ploidy is specified but the length of alleles do not match length of chromosome. val size: " );
#else
							DBG_ASSERT( src.size() == pop.totNumLoci(), ValueError,
								"Ploidy is specified but the length of alleles do not match length of chromosome. val size: " );
#endif
							copy(src.begin(), src.end(), pop.ind(ind).genoBegin(m_atPloidy));
						}
					}
					else						  // with m_loci
					{
						if(m_atPloidy==-1)		  // all copied of chromosome
						{
#ifdef SIMUMPI
							DBG_ASSERT(noMarker || src.size() == m_atLoci.size() ||
								(src.size() == m_atLoci.size() * pop.ploidy() ),
								ValueError, "Length of value does not atLoci size");
#else
							DBG_ASSERT(src.size() == m_atLoci.size() ||
								(src.size() == m_atLoci.size() * pop.ploidy() ),
								ValueError, "Length of value does not atLoci size");
#endif
#ifdef SIMUMPI
							if(noMarker)
								return true;
#endif
							for(size_t loc = 0; loc != srcSz ;++loc)
								*(pop.indGenoBegin(ind) + m_atLoci[loc%lociSz] + loc/lociSz*totNumLoci ) = src[loc];
						}
						else					  // one of the copies.
						{
							DBG_ASSERT( src.size() == m_atLoci.size(), ValueError,
								"Ploidy is specified but the length of alleles do not match length of given allele array.");
#ifdef SIMUMPI
							if(noMarker)
								return true;
#endif
							for(size_t loc = 0; loc != srcSz ;++loc)
								*(pop.ind(ind).genoBegin(m_atPloidy) +
								m_atLoci[loc%lociSz] + loc/lociSz*totNumLoci ) = src[loc];
						}
					}
					if( m_sex.empty())
					{
						if( rng().randUniform01() < m_maleFreq )
							pop.ind(ind).setSex( Male );
						else
							pop.ind(ind).setSex( Female );
					}
					else
					{
						pop.ind(ind).setSex( nextSex() );
					}							  // set sex
				}
			}
		}
		else									  // use proportion.
		{
			for(size_t rg = 0; rg < m_ranges.size(); ++rg)
			{
				ULONG left = m_ranges[rg][0], right = m_ranges[rg][1];

				DBG_FAILIF( left > pop.popSize() || right > pop.popSize() || left > right ,
					ValueError, "Invaid m_ranges boundary: " + toStr(left) + " - " + toStr(right));

				size_t srcSz = m_value[0].size();
				size_t lociSz=m_atLoci.size();
#ifdef SIMUMPI
				size_t totNumLoci = pop.localNumLoci();
#else
				size_t totNumLoci = pop.totNumLoci();
#endif
				Weightedsampler ws(rng(), m_proportion);

				for (ULONG ind = left; ind < right; ++ind)
				{
#ifdef SIMUMPI
												  // whole chromosome or geno
					if(m_atLoci.empty() && !noMarker)
#else
						if(m_atLoci.empty())	  // whole chromosome or geno
#endif
					{
												  // by ploidy
						if( srcSz == totNumLoci)
						{
							if(m_atPloidy==-1)
							{
								for(UINT p = 0; p<pop.ploidy(); ++p)
								{
									UINT idx = ws.get();
									copy(m_value[idx].begin(), m_value[idx].end(),
										pop.indGenoBegin(ind)+p*totNumLoci);
								}
							}
							else
							{
								UINT idx = ws.get();
								copy(m_value[idx].begin(), m_value[idx].end(),
									pop.indGenoBegin(ind)+m_atPloidy*totNumLoci);
							}
						}
						else					  // whole geno
						{
							UINT idx = ws.get();
							if(m_atPloidy==-1)
								copy(m_value[idx].begin(), m_value[idx].end(), pop.indGenoBegin(ind));
							else				  // only one copy of chromosome
								copy(m_value[idx].begin(), m_value[idx].end(),
									pop.ind(ind).genoBegin(m_atPloidy));
						}
					}
					else						  /// atLoci is in effect
					{
#ifdef SIMUMPI						
						if (noMarker)
							return true;
#endif							
						if( srcSz == lociSz )	  // one by one
						{
							if(m_atPloidy==-1)
							{
								for(UINT p = 0; p<pop.ploidy(); ++p)
								{
									UINT idx = ws.get();
									for(size_t loc = 0; loc != srcSz ;++loc)
										*(pop.indGenoBegin(ind) + p*totNumLoci + m_atLoci[loc]) = m_value[idx][loc];
								}
							}
							else
							{
								UINT idx = ws.get();
								for(size_t loc = 0; loc != srcSz ;++loc)
									*(pop.indGenoBegin(ind) + m_atPloidy*totNumLoci +
									m_atLoci[loc]) = m_value[idx][loc];
							}

						}
						else					  // who geno (at loci .. though)
						{
							if(m_atPloidy==-1)
							{
								UINT idx = ws.get();
								for(size_t loc = 0; loc != srcSz ;++loc)
									*(pop.indGenoBegin(ind) + m_atLoci[loc%lociSz] +
									loc/lociSz*totNumLoci ) = m_value[idx][loc];
							}
							else
							{
								UINT idx = ws.get();
								for(size_t loc = 0; loc != srcSz ;++loc)
									*(pop.ind(ind).genoBegin(m_atPloidy)
									+ m_atLoci[loc%lociSz] +
									loc/lociSz*totNumLoci ) = m_value[idx][loc];
							}
						}
					}
					if( m_sex.empty())
					{
						if( rng().randUniform01() < m_maleFreq )
							pop.ind(ind).setSex( Male );
						else
							pop.ind(ind).setSex( Female );
					}
					else
					{
						pop.ind(ind).setSex( nextSex() );
					}							  // set sex
				}
			}

		}
		return true;
	}

	bool pyInit::apply(population& pop)
	{
		this->initSexIter();

		for(UINT al = 0, alEnd=pop.totNumLoci(); al < alEnd; ++al)
		{
#ifdef SIMUMPI
			if (mpiRank() != pop.rankOfLocus(al))
				continue;
#endif
			for(UINT sp=0, numSP=pop.numSubPop(); sp < numSP; ++sp)
			{
				for(ULONG it=0, itEnd=pop.subPopSize(sp); it<itEnd; ++it)
				{
					for(UINT p=0, pEnd=pop.ploidy(); p<pEnd; ++p)
					{
						int resInt;
						PyCallFunc3(m_func, "(iii)", al, p, sp, resInt, PyObj_As_Int);
						pop.ind(it,sp).setAllele( static_cast<Allele>(resInt), al, p);
					}
				}
			}
		}
		// initialize sex
		// call randUnif once for each individual
		// (initialize allele need to call randUnif for each locus
		if(this->m_sex.empty())
		{
			for (population::IndIterator it = pop.indBegin(), itEnd=pop.indEnd();
				it != itEnd; ++it)
			{
				if( rng().randUniform01() < this->m_maleFreq )
					it->setSex( Male );
				else
					it->setSex( Female );
			}
		}
		else
		{
			for (population::IndIterator it = pop.indBegin(), itEnd=pop.indEnd();
				it != itEnd; ++it)
			{
				it->setSex( this->nextSex() );
			}
		}
		return true;
	}

}
