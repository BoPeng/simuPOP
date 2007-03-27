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

#include "stator.h"

namespace simuPOP
{
	bool pyEval::apply(population& pop)
	{
		if(m_exposePop)
		{
			PyObject* popObj = pyPopObj(static_cast<void*>(&pop));
			if( popObj == NULL)
				throw SystemError("Could not expose population pointer. Compiled with the wrong version of SWIG? ");

			// set dictionary variable pop to this object
			pop.setVar("pop", popObj);
		}

		m_expr.setLocalDict(pop.dict());
		string res = m_expr.valueAsString();

		if( ! this->noOutput() )
		{
			ostream & out = this->getOstream(pop.dict());
			out << res;
			this->closeOstream();
		}
		return true;
	}

	string haploKey(const vectori& seq)
	{
		string key = "{'"+toStr(seq[0]);

		for(size_t i=1; i< seq.size(); ++i)
			key += toStr("-") + toStr(seq[i]);

		return key + "'}";
	}

	bool stat::apply(population& pop)
	{
		return (m_popSize.apply(pop) &&
			m_numOfMale.apply(pop) &&
			m_numOfAffected.apply(pop) &&
			m_alleleFreq.apply(pop) &&
			m_heteroFreq.apply(pop) &&
			m_expHetero.apply(pop) &&
			m_genoFreq.apply(pop) &&
			m_haploFreq.apply(pop) &&
			m_LD.apply(pop) &&
			m_association.apply(pop) &&
			m_Fst.apply(pop) &&
			m_relatedness.apply(pop));
	}

	bool statPopSize::apply(population& pop)
	{
		if( !m_isActive)
			return true;

		UINT numSP = pop.numSubPop();
		ULONG popSize = pop.popSize();

		pop.setIntVar(numSubPop_String, numSP);
		pop.setIntVar(popSize_String, popSize);

		// type mismatch, can not use subPopSizes() directly.
		vectori spSize(numSP);
		for(size_t sp=0; sp < numSP; ++sp)
			spSize[sp] = pop.subPopSize(sp);

		pop.setIntVectorVar(subPopSize_String, spSize);

		for( size_t sp=0; sp < numSP; ++sp)
			pop.setIntVar(subPopVar_String(sp, popSize_String), spSize[sp]);
		return true;
	}

	bool statNumOfMale::apply(population& pop)
	{
		if(m_numOfMale.empty())
			return true;

		UINT numSP = pop.numSubPop();
		m_numOfMale.resize(numSP+1);
		m_numOfFemale.resize(numSP+1);

		ULONG numOfMale=0;

		for( size_t sp=0; sp < numSP; ++sp)
		{
			ULONG n = 0;
			for(population::IndIterator it = pop.indBegin(sp), itEnd=pop.indEnd(sp);
				it < itEnd; ++it)
			{
				if(it->sex() == Male)
					n++;
			}
			numOfMale += n;
			m_numOfMale[sp] = n;

			pop.setIntVar(subPopVar_String(sp, numOfMale_String), n);

			n = pop.subPopSize(sp) - n;
			m_numOfFemale[sp] = n;

			pop.setIntVar( subPopVar_String(sp, numOfFemale_String), n);
		}
		pop.setIntVar( numOfMale_String, numOfMale);
		pop.setIntVar( numOfFemale_String, pop.popSize() - numOfMale);
		m_numOfMale[numSP] = numOfMale;
		m_numOfFemale[numSP] = pop.popSize() - numOfMale;
		return true;
	}

	bool statNumOfAffected::apply(population& pop)
	{
		if( m_numOfAffected.empty() )
			return true;

		ULONG numOfAffected=0;
		UINT numSP = pop.numSubPop();
		m_numOfAffected.resize(numSP+1);
		m_numOfUnaffected.resize(numSP+1);

		for( size_t sp=0; sp < numSP; ++sp)
		{
			ULONG n = count_if(pop.indBegin( sp ), pop.indEnd( sp ),
				isAffected<individual>());
			numOfAffected += n;
			m_numOfAffected[sp] = n;
			pop.setIntVar(subPopVar_String(sp, numOfAffected_String), n);
			pop.setDoubleVar( subPopVar_String(sp, propOfAffected_String),
				(double)(n)/pop.subPopSize(sp));

			n = pop.subPopSize(sp) - n;
			m_numOfUnaffected[sp] = n;

			pop.setIntVar(subPopVar_String(sp, numOfUnaffected_String), n);
			pop.setDoubleVar( subPopVar_String(sp, propOfUnaffected_String),
				(double)(n)/pop.subPopSize(sp));
		}
		pop.setIntVar( numOfAffected_String, numOfAffected);
		pop.setIntVar( numOfUnaffected_String, pop.popSize() - numOfAffected);
		pop.setDoubleVar( propOfAffected_String, (double)(numOfAffected)/pop.popSize());
		pop.setDoubleVar( propOfUnaffected_String, (double)(pop.popSize() - numOfAffected)/pop.popSize());
		m_numOfAffected[numSP] = numOfAffected;
		m_numOfUnaffected[numSP] = pop.popSize() - numOfAffected;
		return true;
	}

	bool statAlleleFreq::apply(population& pop)
	{
		if( m_atLoci.empty())
			return true;

		pop.removeVar(NumOfAlleles_String);
		pop.removeVar(AlleleNum_String);
		pop.removeVar(AlleleFreq_String);

		UINT numSP = pop.numSubPop();
		UINT numLoci = m_atLoci.size();

		UINT len = numSP==1?1:(numSP+1);
		// if not initialized or m_atLoci/numSP changes
		if( m_alleleNum.size() != len )
		{
			for( size_t i=0; i < numLoci;  ++i)
				DBG_FAILIF( static_cast<UINT>(m_atLoci[i]) >= pop.totNumLoci(),
					IndexError, "locus index (" + toStr(m_atLoci[i])
					+ ") out of range of 0 - " + toStr(pop.totNumLoci()-1));

			m_alleleNum.resize(len);
			m_alleleFreq.resize(len);

			m_numOfAlleles.resize(len);
			for(size_t i = 0; i < len; ++i)
			{
				m_alleleNum[i].resize(pop.totNumLoci());
				m_alleleFreq[i].resize(pop.totNumLoci());
				m_numOfAlleles[i].resize(pop.totNumLoci());
				fill(m_numOfAlleles[i].begin(), m_numOfAlleles[i].end(), 0);
			}
		}

		string varname;

		for(size_t i = 0; i < numLoci; ++i)
		{
			UINT loc = m_atLoci[i];

			vectori& sum = m_alleleNum.back()[loc];
			fill(sum.begin(), sum.end(), 0);

			// for each subpopulation
			for( UINT sp=0; sp < numSP;  ++sp)
			{
				vectori& num = m_alleleNum[sp][loc];
				// clear all current values
				fill(num.begin(), num.end(),0);

				// go through all alleles
				for(GappedAlleleIterator a=pop.alleleBegin(loc, sp, false),
					aEnd=pop.alleleEnd(loc, sp, false); a != aEnd; ++a)
				{
					if( AlleleUnsigned(*a) >= num.size() )
						num.resize(*a+1, 0);
					num[*a]++;
				}

				// add this number to overall num
				// calculate frequency
				// if there is only one sp, no need to do so.
				if(numSP > 1 )
				{
					if(sum.size() < num.size())
						sum.resize(num.size(), 0);
					for(size_t e=0, eEnd=num.size(); e < eEnd; ++e)
						sum[e] += num[e];
				}

				vectorf& freq = m_alleleFreq[sp][loc];
				freq.resize( num.size(), 0.);
				double dy = pop.subPopSize(sp)*pop.ploidy();
				for(size_t e=0, eEnd=num.size(); e < eEnd; ++e)
					freq[e] = static_cast<double>(num[e])/dy;

				// post result at this locus
				if(m_ifPost[i])
				{
					varname = subPopVar_String(sp, AlleleNum_String) + "[" + toStr(loc) + "]";
					PyObject * d = pop.setIntVectorVar(varname, num);

					// do not need a separate result
					if(numSP == 1)
					{
						varname = toStr(AlleleNum_String) + "[" + toStr(loc) + "]";
						Py_INCREF(d);
						pop.setVar(varname, d);
					}

					varname = subPopVar_String(sp, AlleleFreq_String) + "[" + toStr(loc) + "]";
					d = pop.setDoubleVectorVar(varname, freq);

					// do not need a separate result
					if(numSP == 1)
					{
						varname = toStr(AlleleFreq_String) + "[" + toStr(loc) + "]";
						Py_INCREF(d);
						pop.setVar(varname, d);
					}
				}								  // post

				// set numOfAlleles if necessary
				m_numOfAlleles[sp][loc] = count_if( num.begin(), num.end(),
					bind2nd(std::greater<int>(), 0));
			}									  // subpop

			if(numSP > 1 )						  // calculate sum and post overall result
			{
				// summary?
				vectorf& freq = m_alleleFreq.back()[loc];
				freq.resize( sum.size());
				double dy = pop.popSize() * pop.ploidy();
				for(size_t e=0, eEnd=sum.size(); e < eEnd; ++e)
					freq[e] = static_cast<double>(sum[e])/dy;

				if(m_ifPost[i])
				{
					varname = string(AlleleNum_String) + "[" + toStr(loc) + "]";
					pop.setIntVectorVar(varname, sum);

					varname = string(AlleleFreq_String) + "[" + toStr(loc) + "]";
					pop.setDoubleVectorVar(varname, freq);
				}

				// set numOfAlleles if necessary
				m_numOfAlleles.back()[loc] = count_if( sum.begin(), sum.end(),
					bind2nd(std::greater<int>(), 0));
			}
		}										  // all loci

		// post number of alleles
		if( accumulate(m_ifPost.begin(), m_ifPost.end(), 0) > 0 )
		{
			// post number of alleles
			for(UINT sp = 0; sp < numSP; ++sp)
			{
				PyObject* d = pop.setIntVectorVar( subPopVar_String(sp, NumOfAlleles_String),
					m_numOfAlleles[sp]);
				if(numSP == 1)
				{
					Py_INCREF(d);
					pop.setVar(NumOfAlleles_String, d);
				}
			}
			if(numSP>1)
			{
				pop.setIntVectorVar(NumOfAlleles_String, m_numOfAlleles.back());
			}
		}
		return true;
	}

	bool statHeteroFreq::apply(population& pop)
	{
		if( m_atLoci.empty())
			return true;

		pop.removeVar(HeteroNum_String);
		pop.removeVar(HeteroFreq_String);
		pop.removeVar(AllHeteroNum_String);
		pop.removeVar(AllHeteroFreq_String);
		pop.removeVar(HomoNum_String);
		pop.removeVar(HomoFreq_String);

		UINT numSP = pop.numSubPop();
		UINT numLoci = m_atLoci.size();

		// may be resizing for different replicate of populations.
		// if not initialized or m_atLoci/numSP changes
		if(m_heteroNum.size() != (numSP+1)*numLoci)
		{
			m_heteroNum.resize((numSP+1)*numLoci);
			m_heteroFreq.resize((numSP+1)*numLoci);
			m_homoNum.resize( numSP+1 );
			m_homoFreq.resize( numSP+1 );
		}

		string varname;
		ULONG popSize = pop.popSize();

		for( size_t i = 0; i < numLoci; ++i)
		{
			UINT loc = m_atLoci[i];
			DBG_DO(DBG_STATOR, cout << "Counting heterozygotes at locus " << loc << endl);

			vectori& sum = m_heteroNum[resIdx(i)];
			fill(sum.begin(), sum.end(), 0);
			int sumAll = 0;

			// for each subpopulation
			for( UINT sp=0; sp < numSP;  ++sp)
			{
				vectori& num = m_heteroNum[resIdx(i,sp)];
				fill(num.begin(), num.end(), 0 );
				int numAll = 0;

				// go through all alleles
				//?>> \todo here we assume diploid population
				for(GappedAlleleIterator a=pop.alleleBegin(loc, sp, false),
					aEnd=pop.alleleEnd(loc, sp, false);
					a != aEnd; a+=2)
				{
					if( AlleleUnsigned(*a) >= num.size() )
						num.resize(*a+1);

					if( AlleleUnsigned(*(a+1)) >= num.size() )
						num.resize(*(a+1)+1);

					if( *a != *(a+1) )
					{
						num[*a]++;
						num[*(a+1)]++;
						numAll++;
					}
				}

				sumAll += numAll;
				// add this number to overall num
				// calculate frequency
				if(numSP > 1)
				{
					if(sum.size() < num.size())
						sum.resize(num.size());
					for(size_t e=0, eEnd=num.size(); e < eEnd; ++e)
						sum[e] += num[e];
				}

				vectorf& freq = m_heteroFreq[resIdx(i,sp)];
				freq.resize( num.size());
				for(size_t e=0, eEnd=num.size(); e < eEnd; ++e)
					freq[e] = static_cast<double>(num[e])/pop.subPopSize(sp);

				// set variable.
				if(m_ifPost[i] && m_postHetero)
				{
					varname =  subPopVar_String(sp, HeteroNum_String) + "[" + toStr(loc) + "]";
					PyObject* d = pop.setIntVectorVar(varname, num);
					if(numSP == 1)
					{
						Py_INCREF(d);
						varname =  toStr(HeteroNum_String) + "[" + toStr(loc) + "]";
						pop.setVar(varname, d);
					}

					varname = subPopVar_String(sp, HeteroFreq_String) + "[" + toStr(loc) + "]";
					d = pop.setDoubleVectorVar(varname, freq);

					if(numSP == 1)
					{
						Py_INCREF(d);
						varname =  toStr(HeteroFreq_String) + "[" + toStr(loc) + "]";
						pop.setVar(varname, d);
					}

					// overall hetero
					varname =  subPopVar_String(sp, AllHeteroNum_String) + "[" + toStr(loc) + "]";
					d = pop.setIntVar(varname, numAll);
					if(numSP == 1)
					{
						Py_INCREF(d);
						varname =  toStr(AllHeteroNum_String) + "[" + toStr(loc) + "]";
						pop.setVar(varname, d);
					}

					varname = subPopVar_String(sp, AllHeteroFreq_String) + "[" + toStr(loc) + "]";
					d = pop.setDoubleVar(varname, static_cast<double>(numAll)/pop.subPopSize(sp));

					if(numSP == 1)
					{
						Py_INCREF(d);
						varname =  toStr(AllHeteroFreq_String) + "[" + toStr(loc) + "]";
						pop.setVar(varname, d);
					}
				}
			}

			if(numSP > 1 && m_postHetero)
			{
				vectorf& freq = m_heteroFreq[resIdx(i)];
				freq.resize( sum.size());
				for(size_t e=0, eEnd=sum.size(); e < eEnd; ++e)
					freq[e] = static_cast<double>(sum[e])/popSize;

				if(m_ifPost[i])
				{								  //
					varname = string(HeteroNum_String) + "[" + toStr(loc) + "]";
					pop.setIntVectorVar(varname, sum);

					varname = string(HeteroFreq_String) + "[" + toStr(loc) + "]";
					pop.setDoubleVectorVar(varname, freq);

					varname = string(AllHeteroNum_String) + "[" + toStr(loc) + "]";
					pop.setIntVar(varname, sumAll);

					varname = string(AllHeteroFreq_String) + "[" + toStr(loc) + "]";
					pop.setDoubleVar(varname, static_cast<double>(sumAll)/pop.popSize());
				}
			}									  // whole population
		}										  // for all loci

		if( m_postHomo )
		{
			for( size_t i = 0; i < numLoci; ++i)
			{
				UINT loc = m_atLoci[i];

				if(loc+1 >= m_homoFreq[0].size())
				{
					for( UINT sp=0; sp < numSP+1;  ++sp)
					{
						m_homoFreq[sp].resize(loc+1, 0.0);
						m_homoNum[sp].resize(loc+1, 0);
					}
				}

				// calculate homoNum
				for( UINT sp=0; sp < numSP;  ++sp)
				{
					m_homoNum[sp][loc] = pop.subPopSize(sp) - m_heteroNum[resIdx(i, sp)][0];
					m_homoFreq[sp][loc] = (double)(m_homoNum[sp][loc])/pop.subPopSize(sp);
				}
				m_homoNum[numSP][loc] = pop.popSize() - m_heteroNum[resIdx(i)][0];
				m_homoFreq[numSP][loc] = (double)(m_homoNum[numSP][loc])/pop.popSize();
			}									  // all loci
			// post result
			for( UINT sp=0; sp < numSP; ++sp)
			{
				pop.setIntVectorVar(subPopVar_String(sp,HomoNum_String),
					m_homoNum[sp]);
				pop.setDoubleVectorVar(subPopVar_String(sp,HomoFreq_String),
					m_homoFreq[sp]);
			}
			pop.setIntVectorVar(HomoNum_String, m_homoNum[numSP]);
			pop.setDoubleVectorVar(HomoFreq_String, m_homoFreq[numSP]);
		}
		return true;
	}

	bool statExpHetero::apply(population& pop)
	{
		if( m_atLoci.empty())
			return true;

		pop.removeVar(ExpHetero_String);

		UINT numSP = pop.numSubPop();
		UINT numLoci = m_atLoci.size();

		if(m_expHetero.size() != numSP+1)
			m_expHetero.resize(numSP+1 );

		for( size_t i = 0; i < numLoci; ++i)
		{
			UINT loc = m_atLoci[i];

			if(loc+1 >= m_expHetero[0].size())
			{
				for( UINT sp=0; sp < numSP+1;  ++sp)
					m_expHetero[sp].resize(loc+1, 0.0);
			}

			// for each subpopulation
			for( UINT sp=0; sp < numSP;  ++sp)
			{
				// calculate expected heterozygosity
				// get allele frequency
				vectorf& af = m_alleleFreq.alleleFreqVec(loc, sp);
				double expHeter=1;
				// 1-sum pi^2
				for(int al = 0, alEnd=af.size() ; al < alEnd; al++)
					expHeter -= af[al]*af[al];

				m_expHetero[sp][loc] = expHeter;
			}

			vectorf& af = m_alleleFreq.alleleFreqVec(loc);
			double expHeter=1;
			// 1-sum pi^2
			for(int al = 0, alEnd = af.size(); al < alEnd; al++)
				expHeter -= af[al]*af[al];

			m_expHetero[numSP][loc] = expHeter;
		}

		// post result
		for( UINT sp=0; sp < numSP; ++sp)
			pop.setDoubleVectorVar(subPopVar_String(sp,ExpHetero_String),
				m_expHetero[sp]);
		pop.setDoubleVectorVar(ExpHetero_String,
			m_expHetero[numSP]);
		return true;
	}

	bool statGenoFreq::apply(population& pop)
	{
		if( m_atLoci.empty())
			return true;

		pop.removeVar(GenotypeNum_String);
		pop.removeVar(GenotypeFreq_String);

		UINT numSP = pop.numSubPop();
		ULONG popSize = pop.popSize();

		for( size_t i=0, iEnd = m_atLoci.size(); i < iEnd;  ++i)
		{
			if(static_cast<UINT>(m_atLoci[i]) >= pop.totNumLoci() )
				throw IndexError("Absolute locus index "
					+ toStr( m_atLoci[i]) + " is out of range of 0 ~ "
					+ toStr( pop.totNumLoci() - 1));
		}

		// first remove genoNum that may be set by previous count.
		// for example if genoNum[a][10] was set but this time there
		// is no allele 10, we will not reset genoNum[a][10] ...
		pop.removeVar( GenotypeNum_String);
		pop.removeVar( GenotypeFreq_String);
		for( UINT sp=0; sp < numSP;  ++sp)
		{
			// remove genoNum, genoFreq that may be set by previous run.
			pop.removeVar( subPopVar_String(sp, GenotypeNum_String));
			pop.removeVar( subPopVar_String(sp, GenotypeFreq_String));
		}

		string varname;

		// deal with genotype
		for( size_t i=0, iEnd = m_atLoci.size(); i < iEnd;  ++i)
		{
			// for each locus, we need to use a vector of dictionaries.
			vector<intDict> sum;

			int loc = m_atLoci[i];

#ifndef BINARYALLELE
			Allele a, b;
#else
			unsigned short a, b;
#endif

			// for each subpopulation
			for( UINT sp=0; sp < numSP;  ++sp)
			{
				DBG_DO(DBG_STATOR, cout << "Counting genotypes at locus " <<
					loc << " subPop " << sp << endl);

				vector<intDict> num;

				/// go through a single allele for all individual, all diploid
				for( GappedAlleleIterator it = pop.alleleBegin(loc, sp, false),
					itEnd = pop.alleleEnd(loc, sp, false); it != itEnd;  it+=2 )
				{
					a = *it;
					b = *(it+1);

					if( !m_phase && a > b )
						std::swap(a,b);

					if( a >= num.size() )
						num.resize(a+1);

					num[a][b]++;

					if( a >= sum.size() )
						sum.resize(a+1);

					sum[a][b]++;
				}

				// register values for this subpopulation
				for(a=0; a < num.size(); ++a)
				{
					/// need to replace previous values
					// if( num[a].empty() )
					//   continue;

					// empty dictionary should be allowed
					varname = subPopVar_String(sp, GenotypeNum_String) +
						+ "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
					pop.setIntDictVar( varname, num[a] );

					// apply frequency
					for(intDict::iterator it=num[a].begin(), itEnd=num[a].end(); it!=itEnd; ++it)
						it->second = it->second/ pop.subPopSize(sp);

					varname =  subPopVar_String(sp, GenotypeFreq_String) +
						+ "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
					pop.setIntDictVar( varname, num[a] );
				}
			}

			for(a=0; a<sum.size(); ++a)
			{
				// if( sum[a].empty() )
				//  continue;

				// empty dictionary should be allowed
				varname = toStr(GenotypeNum_String) + "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
				pop.setIntDictVar( varname, sum[a] );

				// apply frequency
				for(intDict::iterator it=sum[a].begin(), itEnd=sum[a].end(); it!=itEnd; ++it)
					it->second = it->second/ popSize;

				varname =  toStr(GenotypeFreq_String) + "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
				pop.setIntDictVar( varname, sum[a] );
			}
		}
		return true;
	}

	bool statHaploFreq::apply(population& pop)
	{
		if( m_haplotypes.empty())
			return true;

		pop.removeVar(HaplotypeNum_String);
		pop.removeVar(HaplotypeFreq_String);

		UINT nHap = m_haplotypes.size();

		// first time?
		if( m_haploNum.size() != (pop.numSubPop()+1)*nHap)
		{
			for( size_t i = 0; i < nHap;  ++i)
			{
				vectori& haplotype = m_haplotypes[i];
				size_t sz = haplotype.size();

				if( sz == 0 )
					throw ValueError("has to specify some haplotype.");
				else if( sz == 1 )
					throw ValueError("Haplotype must contain alleles from more than one locus.");
			}
			m_haploNum.resize( (pop.numSubPop()+1)*nHap);
			m_haploFreq.resize( (pop.numSubPop()+1)*nHap);
		}

		UINT numSP = pop.numSubPop();

		DBG_DO(DBG_STATOR, cout << "Counting haplotypes" << endl);

		// clear all statistics
		for(size_t h = 0; h < nHap*(numSP+1); ++h)
		{
			m_haploNum[h].clear();
			m_haploFreq[h].clear();
		}

		// for each subpopulation
		for( UINT sp=0; sp < numSP;  ++sp)
		{
			for(size_t h = 0; h < nHap; ++h)
			{
				vectori& haplotype = m_haplotypes[h];
				map< vectori, UINT>& count = m_haploNum[h + sp*nHap];
				map< vectori, UINT>& sum = m_haploNum[h + numSP*nHap];

				size_t sz = haplotype.size();

				vectori sampleHap(sz);

				for( GappedAlleleIterator it = pop.alleleBegin(0, sp, false),
					itEnd = pop.alleleEnd(0, sp, false); it != itEnd;  ++it )
				{
					for( size_t hap=0; hap < sz; ++hap)
						sampleHap[hap] = *(it.ptr()+haplotype[hap]);

					// add sampleHap count
					count[sampleHap] ++;
					sum[sampleHap] ++;
				}
			}
		}
		// finish count

		// calculate haploFreq,
		// for each subpopulation
		for( UINT sp=0; sp < numSP;  ++sp)
		{
			// record both num and freq
			string varNumName =  subPopVar_String(sp, HaplotypeNum_String);
			string varFreqName =  subPopVar_String(sp, HaplotypeFreq_String);
			for(size_t h = 0; h < nHap; ++h)
			{
				vectori& haplotype = m_haplotypes[h];
				map< vectori, UINT>& count = m_haploNum[h + sp*nHap];
				map< vectori, double>& freq = m_haploFreq[h + sp*nHap];

				for(map<vectori, UINT>::iterator it = count.begin(); it!=count.end(); ++it)
				{
					freq[ it->first] = double(it->second)/(pop.subPopSize(sp)*pop.ploidy());
					if(m_ifPost[h])
					{
						pop.setIntVar( varNumName + haploKey(haplotype) + haploKey( it->first) , it->second);
						pop.setDoubleVar( varFreqName + haploKey(haplotype) + haploKey( it->first) ,
							double(it->second)/(pop.subPopSize(sp)*pop.ploidy()));
					}
				}
			}
		}
		// whole population
		for(size_t h = 0; h < nHap; ++h)
		{
			vectori& haplotype = m_haplotypes[h];
			map< vectori, UINT>& count = m_haploNum[h + numSP*nHap];
			map< vectori, double>& freq = m_haploFreq[h+numSP*nHap];

			for(map<vectori, UINT>::iterator it = count.begin(); it!=count.end(); ++it)
			{
				freq[ it->first] = double(it->second)/(pop.popSize()*pop.ploidy());
				if(m_ifPost[h])
				{
					pop.setIntVar( HaplotypeNum_String + haploKey(haplotype) + haploKey( it->first) , it->second);
					pop.setDoubleVar( HaplotypeFreq_String + haploKey(haplotype) + haploKey( it->first) ,
						double(it->second)/(pop.popSize()*pop.ploidy()));
				}
			}
		}
		return true;
	}

	// this function calculate single-allele LD measures
	// D, D_p and r2 are used to return calculated values.
	// LD for subpopulation sp is calculated if subPop is true
	void statLD::calculateLD(const vectori & hapLoci, const vectori & hapAlleles, UINT sp, bool subPop,
		double & P_A, double & P_B, double & D, double & D_prime, double & r2, double & delta2)
	{
		if(subPop)
		{
			// get haplotype freq from the m_haploFreq object
			double P_AB = m_haploFreq.haploFreq(hapLoci, sp)[hapAlleles];
			// get allele freq from the m_alleleFreq object
			P_A = m_alleleFreq.alleleFreq(hapAlleles[0], hapLoci[0], sp);
			P_B = m_alleleFreq.alleleFreq(hapAlleles[1], hapLoci[1], sp);

			// calculate LD
			D = P_AB - P_A * P_B;
			// calculate LD'
			double D_max = D > 0 ? std::min(P_A*(1-P_B), (1-P_A)*P_B):std::min(P_A*P_B,(1-P_A)*(1-P_B));
			// fcmp_eq is the float comparison operator, which treat (-1e-10, 1e-10) or so as 0 (platform dependent)
			D_prime = fcmp_eq(D_max, 0.)?0.:D/D_max;
			r2 = (fcmp_eq(P_A,0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1))?0.:D*D/P_A/(1-P_A)/P_B/(1-P_B);
			// calculate delta2
			delta2 = (fcmp_eq(P_A,0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1))?0.:pow((P_AB*((1-P_A)-(P_B-P_AB))-(P_A-P_AB)*(P_B-P_AB)), 2)/(P_A*(1-P_A)*P_B*(1-P_B));
			// if environmental variable SIMUDEBUG is set to DBG_STATOR, or
			// if TurnOnDebug(DBG_STATOR) is called in python, the following will be printed.
			DBG_DO(DBG_STATOR, cout << "LD: subpop " << sp << " : P_AB: " << P_AB
				<< " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max <<
				" LD: " << D << " LD': " << D_prime << " r2: " << r2 << " delta2: " << delta2 << endl);

		}
		else
		{
			// whole population
			// get haplotype freq
			double P_AB = m_haploFreq.haploFreq(hapLoci)[hapAlleles];
			P_A = m_alleleFreq.alleleFreq(hapAlleles[0], hapLoci[0]);
			P_B = m_alleleFreq.alleleFreq(hapAlleles[1], hapLoci[1]);

			// calculate LD
			D = P_AB - P_A * P_B;
			// calculate LD'
			double D_max = D > 0 ? std::min(P_A*(1-P_B), (1-P_A)*P_B):std::min(P_A*P_B,(1-P_A)*(1-P_B));
			D_prime = fcmp_eq(D_max, 0)?0:D/D_max;
			r2 = (fcmp_eq(P_A,0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1))?0:D*D/P_A/(1-P_A)/P_B/(1-P_B);
			delta2 = (fcmp_eq(P_A,0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1))?0.:pow((P_AB*((1-P_A)-(P_B-P_AB))-(P_A-P_AB)*(P_B-P_AB)), 2)/(P_A*(1-P_A)*P_B*(1-P_B));

			DBG_DO(DBG_STATOR, cout << "LD: P_AB: " << P_AB
				<< " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max <<
				" LD: " << D << " LD': " << D_prime << " r2: " << r2 << " delta2: " << delta2 << endl);
		}
	}

	// this function is called by stat::apply(pop). It is called
	// after m_alleleFreq.apply(pop) and m_haploFreq.apply(pop) so
	// allele and haplotype frequencies should be available.
	bool statLD::apply(population& pop)
	{
		if( m_LD.empty())
			return true;

		UINT numSP = pop.numSubPop();
		UINT nLD = m_LD.size();
		// used for delta2 which can only be computed for 2 alleles
		vectori numofalleles = m_alleleFreq.numOfAlleles();
		bool valid_delta2 = false;

		// remove previous values.
		pop.removeVar(LD_String);
		pop.removeVar(LDPRIME_String);
		pop.removeVar(R2_String);
		pop.removeVar(DELTA2_String);
		pop.removeVar(AvgLDPRIME_String);
		pop.removeVar(AvgR2_String);
		pop.removeVar(AvgDELTA2_String);
		// also vars at each subpopulations
		for( UINT sp=0; sp < numSP;  ++sp)
		{
			// subPopVar_String is nothing but subPop[sp]['string']
			pop.removeVar( subPopVar_String(sp, LD_String));
			pop.removeVar( subPopVar_String(sp, LDPRIME_String));
			pop.removeVar( subPopVar_String(sp, R2_String));
			pop.removeVar( subPopVar_String(sp, DELTA2_String));
			pop.removeVar( subPopVar_String(sp, AvgLDPRIME_String));
			pop.removeVar( subPopVar_String(sp, AvgR2_String));
			pop.removeVar( subPopVar_String(sp, AvgDELTA2_String));
		}
		for(size_t i=0; i < nLD; ++i)
		{
			// specifying alleles
			if( m_LD[i].size() == 4)
			{
				vectori hapLoci(2);
				vectori hapAlleles(2);

				hapLoci[0] = m_LD[i][0];
				hapLoci[1] = m_LD[i][1];
				if (numofalleles[hapLoci[0]] <= 2 && numofalleles[hapLoci[1]] <= 2)
					valid_delta2 = true;
				hapAlleles[0] = m_LD[i][2];
				hapAlleles[1] = m_LD[i][3];

				for( UINT sp=0; sp < numSP;  ++sp)
				{
					// get LD values, P_A, P_B is ignored.
					double D = 0;
					double D_prime = 0;
					double r2 = 0;
					double delta2 = 0;
					double P_A = 0;
					double P_B = 0;
					calculateLD(hapLoci, hapAlleles, sp, true, P_A, P_B, D, D_prime, r2, delta2);

					// haploKey returns ['a-b'] from array a,b
					pop.setDoubleVar( subPopVar_String(sp, LD_String) +
						haploKey( hapLoci) + haploKey(hapAlleles), D);
					pop.setDoubleVar( subPopVar_String(sp, LDPRIME_String) +
						haploKey( hapLoci) + haploKey(hapAlleles), D_prime);
					pop.setDoubleVar( subPopVar_String(sp, R2_String) +
						haploKey( hapLoci) + haploKey(hapAlleles), r2);
					if (valid_delta2)
						pop.setDoubleVar( subPopVar_String(sp, DELTA2_String) +
							haploKey( hapLoci) + haploKey(hapAlleles), delta2);
					// if numSP == 1, use values for the only subpop as whole population
					if (numSP == 1)
					{
						pop.setDoubleVar(LD_String + haploKey( hapLoci) + haploKey(hapAlleles), D);
						pop.setDoubleVar(LDPRIME_String + haploKey( hapLoci) + haploKey(hapAlleles), D_prime);
						pop.setDoubleVar(R2_String + haploKey( hapLoci) + haploKey(hapAlleles), r2);
						if (valid_delta2)
							pop.setDoubleVar(DELTA2_String + haploKey( hapLoci) + haploKey(hapAlleles), delta2);
					}
				}

				if (numSP > 1)
				{
					// whole population, P_A, P_B is ignored
					double D = 0;
					double D_prime = 0;
					double r2 = 0;
					double delta2 = 0;
					double P_A = 0;
					double P_B = 0;
					calculateLD(hapLoci, hapAlleles, 0, false, P_A, P_B, D, D_prime, r2, delta2);

					pop.setDoubleVar(LD_String + haploKey( hapLoci) + haploKey(hapAlleles), D);
					pop.setDoubleVar(LDPRIME_String + haploKey( hapLoci) + haploKey(hapAlleles), D_prime);
					pop.setDoubleVar(R2_String + haploKey( hapLoci) + haploKey(hapAlleles), r2);
					if (valid_delta2)
						pop.setDoubleVar(DELTA2_String + haploKey( hapLoci) + haploKey(hapAlleles), delta2);
				}
			}
			else
			{
				// No alleles specified, average over all available alleles
				vectori hapLoci(2);
				vectori hapAlleles(2);

				hapLoci[0] = m_LD[i][0];
				hapLoci[1] = m_LD[i][1];
				if (numofalleles[hapLoci[0]] <= 2 && numofalleles[hapLoci[1]] <= 2)
					valid_delta2 = true;
				string hapLociStr = '[' + toStr(hapLoci[0]) + "][" +
					toStr(hapLoci[1]) + ']';

				// find out all alleles
				vectori A_alleles = m_alleleFreq.alleles(hapLoci[0]);
				vectori B_alleles = m_alleleFreq.alleles(hapLoci[1]);

				for( UINT sp=0; sp < numSP;  ++sp)
				{
					double D=0.0, D_prime = 0.0, r2 = 0.0, delta2 = 0.0;
					// iterate through all alleles at locus A and B
					for(vectori::iterator A_ale = A_alleles.begin();
						A_ale != A_alleles.end(); ++A_ale)
					{
						for(vectori::iterator B_ale = B_alleles.begin();
							B_ale != B_alleles.end(); ++B_ale)
						{
							// this is now single allele ...
							hapAlleles[0] = *A_ale;
							hapAlleles[1] = *B_ale;
							double D_ = 0;
							double D_prime_ = 0;
							double r2_ = 0;
							double delta2_ = 0;
							double P_A = 0;
							double P_B = 0;
							calculateLD(hapLoci, hapAlleles, sp, true, P_A, P_B, D_, D_prime_, r2_, delta2_);

							// store allele-specific LD values as well.
							if( m_midValues)
							{
								pop.setDoubleVar( subPopVar_String(sp, LD_String) + haploKey(hapLoci) +
									haploKey(hapAlleles), D_);
								pop.setDoubleVar( subPopVar_String(sp, LDPRIME_String) + haploKey(hapLoci) +
									haploKey(hapAlleles), D_prime_);
								pop.setDoubleVar( subPopVar_String(sp, R2_String) + haploKey(hapLoci) +
									haploKey(hapAlleles), r2_);
								if (valid_delta2)
									pop.setDoubleVar( subPopVar_String(sp, DELTA2_String) + haploKey(hapLoci) +
										haploKey(hapAlleles), delta2_);
								// if only one subpopulation, set the same value for whole population
								if( numSP == 1)
								{
									pop.setDoubleVar( LD_String + haploKey(hapLoci) +
										haploKey(hapAlleles), D_);
									pop.setDoubleVar( LDPRIME_String + haploKey(hapLoci) +
										haploKey(hapAlleles), D_prime_);
									pop.setDoubleVar( R2_String + haploKey(hapLoci) +
										haploKey(hapAlleles), r2_);
									if (valid_delta2)
										pop.setDoubleVar( DELTA2_String + haploKey(hapLoci) +
											haploKey(hapAlleles), delta2_);
								}
							}

							D += P_A*P_B*fabs(D_);
							D_prime += P_A*P_B*fabs(D_prime_);
							r2 += P_A*P_B*r2_;
							if (valid_delta2)
								delta2 += P_A*P_B*delta2_;
							DBG_DO(DBG_STATOR, cout << "Sum D " << D << " D' " << D_prime <<" r2 " << r2 << "delta2" << delta2 << endl);
						}
					}

					// average take average.
					pop.setDoubleVar( subPopVar_String(sp, AvgLD_String) + hapLociStr, D);
					pop.setDoubleVar( subPopVar_String(sp, AvgLDPRIME_String) + hapLociStr, D_prime);
					pop.setDoubleVar( subPopVar_String(sp, AvgR2_String) + hapLociStr, r2);
					if (valid_delta2)
						pop.setDoubleVar( subPopVar_String(sp, AvgDELTA2_String) + hapLociStr, delta2);
					// if numSP == 1, use the single subpop value as whole pop
					if( numSP == 1)
					{
						pop.setDoubleVar( AvgLD_String + hapLociStr, D);
						pop.setDoubleVar( AvgLDPRIME_String + hapLociStr, D_prime);
						pop.setDoubleVar( AvgR2_String + hapLociStr, r2);
						if (valid_delta2)
							pop.setDoubleVar( AvgDELTA2_String + hapLociStr, delta2);
					}
				}

				if(numSP > 1 )
				{
					double D = 0.0, D_prime = 0.0, r2 = 0.0, delta2 = 0.0;
					for(vectori::iterator A_ale = A_alleles.begin();
						A_ale != A_alleles.end(); ++A_ale)
					{
						for(vectori::iterator B_ale = B_alleles.begin();
							B_ale != B_alleles.end(); ++B_ale)
						{
							hapAlleles[0] = *A_ale;
							hapAlleles[1] = *B_ale;
							double D_ = 0;
							double D_prime_ = 0;
							double r2_ = 0;
							double delta2_ = 0;
							double P_A = 0;
							double P_B = 0;
							calculateLD(hapLoci, hapAlleles, 0, false, P_A, P_B, D_, D_prime_, r2_, delta2_);

							if( m_midValues)
							{
								pop.setDoubleVar(LD_String + haploKey(hapLoci) + haploKey(hapAlleles), D_);
								pop.setDoubleVar(LDPRIME_String + haploKey(hapLoci) + haploKey(hapAlleles), D_prime_);
								pop.setDoubleVar(R2_String + haploKey(hapLoci) + haploKey(hapAlleles), r2_);
								if (valid_delta2)
									pop.setDoubleVar(DELTA2_String + haploKey(hapLoci) + haploKey(hapAlleles), delta2_);
							}

							D += P_A * P_B * fabs(D_);
							D_prime += P_A * P_B * fabs(D_prime_);
							r2 += P_A * P_B * r2_;
							if (valid_delta2)
								delta2 += P_A * P_B * delta2_;
							DBG_DO(DBG_STATOR, cout << "Sum D " << D << " D' " << D_prime << " r2 " << r2 << endl);
						}						  // all haplotypes
					}
					pop.setDoubleVar( AvgLD_String + hapLociStr, D);
					pop.setDoubleVar( AvgLDPRIME_String + hapLociStr, D_prime);
					pop.setDoubleVar( AvgR2_String + hapLociStr, r2);
					if (valid_delta2)
						pop.setDoubleVar( AvgDELTA2_String + hapLociStr, delta2);
				}								  // length 2
			}
		}										  // for all LD
		return true;
	}

	bool statAssociation::apply(population& pop)
	{
		if( m_association.empty())
			return true;

		UINT numSP = pop.numSubPop();
		UINT nAssociation = m_association.size();

		// remove previous values.
		pop.removeVar(ChiSq_String);
		pop.removeVar(UCU_String);
		pop.removeVar(CramerV_String);
		// also vars at each subpopulations
		for( UINT sp=0; sp < numSP;  ++sp)
		{
			// subPopVar_String is nothing but subPop[sp]['string']
			pop.removeVar( subPopVar_String(sp, ChiSq_String));
			pop.removeVar( subPopVar_String(sp, UCU_String));
			pop.removeVar( subPopVar_String(sp, CramerV_String));
		}
		for(size_t i=0; i < nAssociation; ++i)
		{
			//
			vectori & hapLoci = m_association[i];

			// find out all alleles
			vectori A_alleles = m_alleleFreq.alleles(hapLoci[0]);
			vectori B_alleles = m_alleleFreq.alleles(hapLoci[1]);

			for( UINT sp=0; sp < numSP;  ++sp)
			{
				UINT as = A_alleles.size();
				UINT bs = B_alleles.size();
				vector<vectorf> cont_table(as+1);
				for(size_t i=0; i <= as; ++i)
					cont_table[i].resize(bs+1);
				double n = static_cast<double>(pop.subPopSize(sp));
				double ChiSq = 0.0, UC_U = 0.0, CramerV = 0.0;
				// initialize last line/column
				for(size_t i=0; i < as; ++i)
					cont_table[i][bs]=0;
				for(size_t j=0; j <= bs; ++j)
					cont_table[as][j]=0;
				// get P_ij
				for(size_t i=0; i<as; ++i)
					for(size_t j=0; j < bs; ++j)
				{
					vectori hapAlleles(2);
					hapAlleles[0] = A_alleles[i];
					hapAlleles[1] = B_alleles[j];
					cont_table[i][j] = m_haploFreq.haploFreq(hapLoci)[hapAlleles];
					cont_table[i][bs] += cont_table[i][j];
					cont_table[as][j] += cont_table[i][j];
					cont_table[as][bs] += cont_table[i][j];
				}
				DBG_ASSERT(fcmp_eq(cont_table[as][bs], 1.), ValueError,
					"Sum of haplotype frequencies is not 1. Association will not be computed.");
				DBG_DO(DBG_STATOR, for(size_t i=0; i <= as; ++i) cout << cont_table[i] << endl);
				// calculate statistics
				for(size_t i=0; i<as; ++i)
					for(size_t j=0; j < bs; ++j)
						ChiSq += pow((n*cont_table[i][j] - n*cont_table[i][bs]*cont_table[as][j]), 2)/n*cont_table[i][bs]*cont_table[as][j];
				DBG_DO(DBG_STATOR, cout << "Chisq " << ChiSq << " sp " << sp << endl);
			}
			if(numSP > 1 )
			{
				UINT as = A_alleles.size();
				UINT bs = B_alleles.size();
				vector<vectorf> cont_table(as+1);
				for(size_t i=0; i <= as; ++i)
					cont_table[i].resize(bs+1);
				double n = static_cast<double>(pop.popSize());
				double ChiSq = 0.0, UC_U =0.0, CramerV = 0.0;
				// get p_ij
				for(size_t i=0; i < as; ++i)
					cont_table[i][bs]=0;
				for(size_t j=0; j <= bs; ++j)
					cont_table[as][j]=0;
				for(size_t i=0; i < as; ++i)
					for(size_t j=0; j < bs; ++j)
				{
					vectori hapAlleles(2);
					hapAlleles[0] = A_alleles[i];
					hapAlleles[1] = B_alleles[j];
					cont_table[i][j] = m_haploFreq.haploFreq(hapLoci)[hapAlleles];
					cont_table[i][bs] += cont_table[i][j];
					cont_table[as][j] += cont_table[i][j];
					cont_table[as][bs] += cont_table[i][j];
				}
				DBG_ASSERT(fcmp_eq(cont_table[as][bs], 1.), ValueError,
					"Sum of haplotype frequencies is not 1. Association will not be computed.");
				DBG_DO(DBG_STATOR, for(size_t i=0; i <= as; ++i) cout << cont_table[i] << endl);
				// calculate statistics
				for(size_t i=0; i<as; ++i)
					for(size_t j=0; j < bs; ++j)
						ChiSq += pow((n*cont_table[i][j] - n*cont_table[i][bs]*cont_table[as][j]), 2)/n*cont_table[i][bs]*cont_table[as][j];
				DBG_DO(DBG_STATOR, cout << "Chisq " << ChiSq << endl);
			}
		}
		return true;
	}

	bool statFst::apply(population& pop)
	{
		if( m_atLoci.empty())
			return true;

		pop.removeVar(Fst_String);
		pop.removeVar(Fis_String);
		pop.removeVar(Fit_String);

		m_Fst.clear();
		m_Fit.clear();
		m_Fis.clear();

		// dicitonary to save values.
		UINT numSP = pop.numSubPop();
		ULONG popSize = pop.popSize();

		// do not save these values now
		double aa = 0., bb = 0., cc = 0.;

		// vector to store p[i]
		vectorf p_i = vectorf(numSP);

		// calculate Fst for each locus
		for(size_t st = 0; st < m_atLoci.size(); ++st)
		{
			int loc = m_atLoci[st];

			DBG_ASSERT( static_cast<size_t>(loc) < pop.totNumLoci(), IndexError,
				"Index out of range of 0 ~ " + toStr(pop.totNumLoci()-1));

			// get all available alleles
			vectori alleles = m_alleleFreq.alleles(loc);

			DBG_DO(DBG_STATOR, cout << "Using alleles " << alleles << endl);

			// n_bar
			double r = numSP;
			double n = popSize;
			double n_bar = n/r;
			vectorlu n_i = pop.subPopSizes();

			// n_c
			double n_c = n;
			for( int i = 0; i < r; ++i)
				n_c  -= n_i[i]*n_i[i]/n;
			n_c /= (r - 1);

			double a=0.0, b=0.0, c=0.0;

			for( vectori::iterator ale = alleles.begin(); ale != alleles.end(); ++ale)
			{
				// p_i
				for(int sp=0; sp < r; ++sp)
					p_i[sp] = m_alleleFreq.alleleFreq(*ale, loc, sp);

				// p_bar
				double p_bar = 0;
				for(int sp = 0; sp < r; ++sp)
					p_bar += n_i[sp] * p_i[sp];
				p_bar /= n;

				// s^2
				double s_2 = 0;
				for(int sp=0; sp < r; ++sp)
					s_2 += n_i[sp]*(p_i[sp] - p_bar)*(p_i[sp]-p_bar);
				s_2 /= (r-1)*n_bar;

				// h_bar
				double h_bar = 0;
				for(int sp=0; sp < r; ++sp)
					h_bar += m_heteroFreq.heteroFreq(*ale, loc, sp) * n_i[sp];
				h_bar /= n;

				// a, b, c
				a += n_bar/n_c *( s_2 - ( p_bar*(1-p_bar) - (r-1.)/r*s_2 - h_bar/4.)/(n_bar -1. ) );
				b += n_bar / (n_bar -1 )* ( p_bar *(1-p_bar) - (r-1)/r * s_2 - (2 * n_bar -1 )/(4.*n_bar)* h_bar );
				c += h_bar /2.;

				DBG_DO(DBG_STATOR, cout << "allele " << *ale << "\tn_c: " << n_c
					<< "\tp_i: " << p_i << "\tp_bar: " << p_bar << "\ts^2: " << s_2 << "\th_bar:"
					<< h_bar << "\ta: " << a << "\tb: " << b << "\tc: " << c << endl);
			}									  // each allele

			DBG_DO(DBG_STATOR, cout << "Fst= " << a/(a+b+c) << endl);

			if( static_cast<size_t>(loc) >= m_Fst.size())
			{
				m_Fst.resize(loc+1,0.);
				m_Fit.resize(loc+1,0.);
				m_Fis.resize(loc+1,0.);
			}
			m_Fst[loc] = fcmp_eq(a+b+c,0.)?0.:(a / ( a+b+c));
			m_Fit[loc] = fcmp_eq(a+b+c,0.)?1.:(1 - c / (a+b+c));
			m_Fis[loc] = fcmp_eq(b+c,0.)?1.:(1 - c /( b+c));

			aa += a;
			bb += b;
			cc += c;
		}
		m_avgFst = fcmp_eq(aa+bb+cc,0.)?0:(aa / ( aa+bb+cc));
		m_avgFit = fcmp_eq(aa+bb+cc,0.)?1.:(1 - cc / (aa+bb+cc));
		m_avgFis = fcmp_eq(aa+bb+cc,0)?1.:(1 - cc /( bb+cc));

		// post results
		pop.setDoubleVectorVar(Fst_String, m_Fst);
		pop.setDoubleVectorVar(Fit_String, m_Fit);
		pop.setDoubleVectorVar(Fis_String, m_Fis);
		pop.setDoubleVar(AvgFst_String, m_avgFst);
		pop.setDoubleVar(AvgFit_String, m_avgFit);
		pop.setDoubleVar(AvgFis_String, m_avgFis);
		return true;
	}

	// relatedness between individuals
	statRelatedness::fraction statRelatedness::relQueller(individual ind1,
		individual ind2)
	{
		matrix& af = m_alleleFreq.alleleFreqAll();
		statRelatedness::fraction res(0.,0.);
		int numScored = 0;
		for(vectori::iterator locus=m_atLoci.begin();
			locus != m_atLoci.end(); ++locus)
		{
			Allele a = ind1.allele(*locus, 0);
			Allele b = ind1.allele(*locus, 1);
			Allele c = ind2.allele(*locus, 0);
			Allele d = ind2.allele(*locus, 1);

			if( a==0 || b==0 || c==0 || d==0 ) continue;

			double s = af[*locus][a] + af[*locus][b] + af[*locus][c] + af[*locus][d];
			double r = (a==c) + (b==c) + (a==d) + (b==d);
			double h = 2 + (a==b) + (c==d);

			res.first += (r - s);
			res.second += (h - s);
			numScored ++;
		}
		// chromosome
		// cout << numScored <<  " " << res.value() << endl;
		if( res.second == 0 || numScored <= m_minScored)
			cout << "Warning: Not enough alleles available to calculate relatedness. Num scored = "
				<< numScored << endl;

		return res;
	}

	statRelatedness::fraction statRelatedness::relLynch(individual ind1,
		individual ind2)
	{
		matrix& af = m_alleleFreq.alleleFreqAll();

		statRelatedness::fraction res(0.,0.);
		int numScored = 0;

		double rel_xy = 0.0;
		double rel_yx = 0.0;
		double weight_xy = 0.0;
		double weight_yx = 0.0;
		for(vectori::iterator locus=m_atLoci.begin();
			locus != m_atLoci.end(); ++locus)
		{
			Allele a = ind1.allele(*locus, 0);
			Allele b = ind1.allele(*locus, 1);
			Allele c = ind2.allele(*locus, 0);
			Allele d = ind2.allele(*locus, 1);

			if( a==0 || b==0 || c==0 || d==0 ) continue;

			double pa = af[*locus][a];
			double pb = af[*locus][b];
			double pc = af[*locus][c];
			double pd = af[*locus][d];

			if( pa < 1e-8 || pb < 1e-8 || pc < 1e-8 || pd < 1e-8) continue;

			double r_xy = ( pa*((b==c)+(b==d)) + pb*((a==c)+(a==d)) - 4*pa*pb ) /
				( (1+(a==b))*(pa+pb) - 4*pa*pb);
			double r_yx = ( pc*((d==a)+(d==b)) + pd*((c==a)+(c==b)) - 4*pc*pd ) /
				( (1+(c==d))*(pc+pd) - 4*pc*pd);
			double w_xy = ( (1+(a==b))*(pa+pb) - 4*pa*pb) / (2*pa*pb);
			double w_yx = ( (1+(c==d))*(pc+pd) - 4*pc*pd) / (2*pc*pd);

			rel_xy += w_xy*r_xy;
			rel_yx += w_yx*r_yx;
			weight_xy += w_xy;
			weight_yx += w_yx;

			numScored ++;
		}										  // all loci
		DBG_FAILIF( numScored <= m_minScored || weight_xy == 0.0 || weight_yx == 0.0,
			ValueError, "Not enough allels to calculated relatedness.");
		res.first = (rel_xy/weight_xy + rel_yx/weight_yx)/2.;
		res.second = 1.;
		return res;
	}

	// IR measure for individual ind at specified locus
	statRelatedness::fraction statRelatedness::relIR(individual ind1, int locus)
	{
		matrix& af = m_alleleFreq.alleleFreqAll();
		statRelatedness::fraction res(0.,0.);
		Allele a = ind1.allele(locus, 0);
		Allele b = ind1.allele(locus, 1);
		double pa = af[locus][a];
		double pb = af[locus][b];

		res.first = 2*(a==b) - pa - pb;
		res.second = 2 - pa - pb;

		if( fcmp_eq( res.second, 0) )
			cout << "Warning: IR value: pa = pb =1. NA will be returned." << endl;
		return res;
	}

	// D2 measure for individual ind at specified locus
	statRelatedness::fraction statRelatedness::relD2(individual ind1, int locus)
	{
		statRelatedness::fraction res(0.,1.);
		Allele a = ind1.allele(locus, 0);
		Allele b = ind1.allele(locus, 1);

		UINT mx = m_alleleFreq.numOfAlleles()[locus];
		res.first = ((a - b)/(mx-2.0))*((a - b)/(mx-2.0));

		return res;
	}

	// REL measure for individual ind at specified locus
	statRelatedness::fraction statRelatedness::relRel(individual ind1,
		individual ind2,  int locus)
	{
		matrix& af = m_alleleFreq.alleleFreqAll();

		statRelatedness::fraction res(0.,0.);
		Allele a = ind1.allele(locus, 0);
		Allele b = ind1.allele(locus, 1);
		Allele c = ind2.allele(locus, 0);
		Allele d = ind2.allele(locus, 1);

		if( a==0 || b==0 || c==0 || d==0 )
			return res;

		double s = af[locus][a] + af[locus][b] + af[locus][c] + af[locus][d];
		double r = (a==c) + (b==c) + (a==d) + (b==d);
		double h = 2 + (a==b) + (c==d);

		res.first += (r - s);
		res.second += (h - s);

		return res;
	}

	// between group i and j if method=REL_Queller and REL_Lynch
	/// for group i and locus j otherwise
	double statRelatedness::groupRelatedness(population& pop, int i, int j, int method)
	{
		statRelatedness::fraction res(0., 0.);

		if( m_useSubPop)
		{
			UINT sp1, sp2;
			if( m_groups[0].empty())
			{
				sp1 = i;
				sp2 = j;
			}
			else
			{
				sp1 = m_groups[0][i];
				sp2 = m_groups[0][j];
			}

			switch(method)
			{
				case REL_Queller:
					// from subpop i and j
					for(population::IndIterator ind1 = pop.indBegin(sp1);
						ind1 != pop.indEnd(sp1); ++ind1)
					{
						for(population::IndIterator ind2 = pop.indBegin(sp2);
							ind2 != pop.indEnd(sp2); ++ind2)
						{
							fraction tmp= relQueller(*ind1, *ind2);
							res.first += tmp.first/tmp.second;
							res.second += 1.;
						}
					}
					return res.first/res.second;
				case REL_Lynch:
					// from subpop i and j
					for(population::IndIterator ind1 = pop.indBegin(sp1);
						ind1 != pop.indEnd(sp1); ++ind1)
					{
						for(population::IndIterator ind2 = pop.indBegin(sp2);
							ind2 != pop.indEnd(sp2); ++ind2)
						{
							fraction tmp= relLynch(*ind1, *ind2);
							res.first += tmp.first/tmp.second;
							res.second += 1.;
						}
					}							  // lynch
					return res.first/res.second;
				case REL_IR:
					for(population::IndIterator ind1 = pop.indBegin(sp1);
						ind1 != pop.indEnd(sp1); ++ind1)
					{
						fraction tmp= relIR(*ind1, j);
						res.first += tmp.first;
						res.second += tmp.second;
					}							  // lynch
					return res.first/res.second;
				case REL_D2:
					for(population::IndIterator ind1 = pop.indBegin(sp1);
						ind1 != pop.indEnd(sp1); ++ind1)
					{
						fraction tmp= relD2(*ind1, j);
						res.first += tmp.first/tmp.second;
						res.second += 1.;
					}							  // lynch
					return res.first/res.second;
				case REL_Rel:
					for(population::IndIterator ind1 = pop.indBegin(sp1);
						ind1 != pop.indEnd(sp1); ++ind1)
					{
						for(population::IndIterator ind2 = ind1 + 1;
							ind2 != pop.indEnd(sp1); ++ind2)
						{
							fraction tmp= relRel(*ind1, *ind2, j);
							res.first += tmp.first;
							res.second += tmp.second;
						}
					}							  // lynch
					return res.first/res.second;
			}									  // switch
		}										  // m_useSubPop
		else
		{
			switch(method)
			{
				case REL_Queller:
					// from specified group i and j
					for(vectori::iterator ind1 = m_groups[i].begin();
						ind1 != m_groups[i].end(); ++ind1)
					{
						for(vectori::iterator ind2 = m_groups[j].begin();
							ind2 != m_groups[j].end(); ++ind2)
						{
							fraction tmp = relQueller(pop.ind(*ind1), pop.ind(*ind2));
							res.first += tmp.first/tmp.second;
							res.second += 1.;
						}
					}
					return res.first/res.second;
				case REL_Lynch:
					for(vectori::iterator ind1 = m_groups[i].begin();
						ind1 != m_groups[i].end(); ++ind1)
					{
						for(vectori::iterator ind2 = m_groups[j].begin();
							ind2 != m_groups[j].end(); ++ind2)
						{
							fraction tmp = relLynch(pop.ind(*ind1), pop.ind(*ind2));
							res.first += tmp.first/tmp.second;
							res.second += 1.;
						}
					}
					return res.first/res.second;
				case REL_IR:
					for(vectori::iterator ind1 = m_groups[i].begin();
						ind1 != m_groups[i].end(); ++ind1)
					{
						fraction tmp = relIR(pop.ind(*ind1), j);
						res.first += tmp.first;
						res.second += tmp.second;
					}
					return res.first/res.second;
				case REL_D2:
					for(vectori::iterator ind1 = m_groups[i].begin();
						ind1 != m_groups[i].end(); ++ind1)
					{
						fraction tmp = relD2(pop.ind(*ind1), j);
						res.first += tmp.first/tmp.second;
						res.second += 1.;
					}
					return res.first/res.second;
				case REL_Rel:
					for(vectori::iterator ind1 = m_groups[i].begin();
						ind1 != m_groups[i].end(); ++ind1)
					{
						for(vectori::iterator ind2 = ind1 + 1;
							ind2 != m_groups[i].end(); ++ind2)
						{
							fraction tmp = relRel(pop.ind(*ind1),
								pop.ind(*ind2), j);
							res.first += tmp.first;
							res.second += tmp.second;
						}
					}
					return res.first/res.second;
			}									  // switch
		}										  // m_useSubPop
		// will never reach here.
		return 0;
	}

	bool statRelatedness::apply(population& pop)
	{
		if(m_groups.empty())
			return true;

		pop.removeVar(Rel_Queller_String);
		pop.removeVar(Rel_Lynch_String);
		pop.removeVar(Rel_IR_String);
		pop.removeVar(Rel_D2_String);
		pop.removeVar(Rel_Rel_String);

		// calculate relatendness values between two groups
		UINT nGroups;
		if(m_useSubPop)
		{
			if( m_groups[0].empty() )
				nGroups = pop.numSubPop();
			else
				nGroups = m_groups[0].size();
		}
		else
			nGroups = m_groups.size();

		UINT nLoci = m_atLoci.size();

		for(size_t m=0; m < m_method.size(); ++m)
		{
			switch(m_method[m])
			{
				case REL_Queller:
					m_relQueller.resize(nGroups);
					for(UINT i=0; i<nGroups; ++i)
						m_relQueller[i].resize(nGroups);

					for(UINT i=0; i<nGroups; ++i)
						for(UINT j=i+1; j<nGroups; ++j)
					{
						m_relQueller[i][j] = groupRelatedness(pop, i, j, m_method[m]);
						m_relQueller[j][i] = m_relQueller[i][j];
					}

					for(UINT i=0; i<nGroups; ++i)
						pop.setDoubleVectorVar(Rel_Queller_String + toStr("[")
							+ toStr(i) + "]", m_relQueller[i]);
					break;
				case REL_Lynch:
					m_relLynch.resize(nGroups);
					for(UINT i=0; i<nGroups; ++i)
						m_relLynch[i].resize(nGroups);

					for(UINT i=0; i<nGroups; ++i)
						for(UINT j=i+1; j<nGroups; ++j)
					{
						m_relLynch[i][j] = groupRelatedness(pop, i, j, m_method[m]);
						m_relLynch[j][i] = m_relLynch[i][j];
					}

					for(UINT i=0; i<nGroups; ++i)
						pop.setDoubleVectorVar(Rel_Lynch_String + toStr("[")
							+ toStr(i) + "]", m_relLynch[i]);
					break;
				case REL_IR:
					m_relIR.resize(nGroups);
					for(UINT i=0; i<nGroups; ++i)
						m_relIR[i].resize(nLoci);

					for(UINT i=0; i<nGroups; ++i)
						for(UINT j=0; j<nLoci; ++j)
							m_relIR[i][j] = groupRelatedness(pop, i, m_atLoci[j], m_method[m]);

					for(UINT i=0; i<nGroups; ++i)
						pop.setDoubleVectorVar(Rel_IR_String + toStr("[")
							+ toStr(i) + "]", m_relIR[i]);
					break;
				case REL_D2:
					m_relD2.resize(nGroups);
					for(UINT i=0; i<nGroups; ++i)
						m_relD2[i].resize(nLoci);

					for(UINT i=0; i<nGroups; ++i)
						for(UINT j=0; j<nLoci; ++j)
							m_relD2[i][j] = groupRelatedness(pop, i, m_atLoci[j], m_method[m]);

					for(UINT i=0; i<nGroups; ++i)
						pop.setDoubleVectorVar(Rel_D2_String + toStr("[")
							+ toStr(i) + "]", m_relD2[i]);
					break;
				case REL_Rel:
					m_relRel.resize(nGroups);
					for(UINT i=0; i<nGroups; ++i)
						m_relRel[i].resize(nLoci);

					for(UINT i=0; i<nGroups; ++i)
						for(UINT j=0; j<nLoci; ++j)
							m_relRel[i][j] = groupRelatedness(pop, i, m_atLoci[j], m_method[m]);

					for(UINT i=0; i<nGroups; ++i)
						pop.setDoubleVectorVar(Rel_Rel_String + toStr("[")
							+ toStr(i) + "]", m_relRel[i]);
					break;
			}
		}										  // for all method
		return true;
	}

}
