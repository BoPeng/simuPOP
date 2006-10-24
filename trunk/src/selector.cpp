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

#include "selector.h"

namespace simuPOP
{
	double mapSelector::indFitness(individual * ind)
	{
		string key;

		for(vectoru::iterator loc=m_loci.begin(); loc!=m_loci.end(); ++loc)
		{
			/// get genotype of ind
			Allele a = ind->allele(*loc, 0);
			Allele b = ind->allele(*loc, 1);

			if( loc != m_loci.begin() )
				key += '|';
			if( ! m_phase && a > b )			  // ab=ba
				key +=  toStr(static_cast<int>(b)) + "-" + toStr(static_cast<int>(a));
			else
				key +=  toStr(static_cast<int>(a)) + "-" + toStr(static_cast<int>(b));
		}

		strDict::iterator pos = m_dict.find(key);

		DBG_ASSERT( pos != m_dict.end(), ValueError,
			"No fitness value for genotype " + key);

		return( pos->second);
	}

	/// currently assuming diploid
	double maSelector::indFitness(individual * ind)
	{
		UINT index = 0;
		bool singleST = m_wildtype.size() == 1;
		for(vectoru::iterator loc=m_loci.begin(); loc!=m_loci.end(); ++loc)
		{
			/// get genotype of ind
			Allele a = ind->allele(*loc, 0);
			Allele b = ind->allele(*loc, 1);

			int numWildtype=0;

			// count number of wildtype
			// this improve the performance a little bit
			if(singleST)
			{
				numWildtype = (a==m_wildtype[0]) + (b==m_wildtype[0]);
			}
			else
			{
				if(find(m_wildtype.begin(), m_wildtype.end(), a) != m_wildtype.end() )
					numWildtype ++;

				if(find(m_wildtype.begin(), m_wildtype.end(), b) != m_wildtype.end() )
					numWildtype ++;
			}

			index = index*3 + 2-numWildtype;
		}
		return m_fitness[index];
	}

	double mlSelector::indFitness(individual * ind)
	{
		if(m_mode == SEL_Multiplicative )
		{
			double fit = 1;
			for(vectorop::iterator s = m_selectors.begin(), sEnd=m_selectors.end();
				s != sEnd; ++s)
			fit *= static_cast<selector* >(*s)->indFitness(ind);
			return fit;
		}
		else if(m_mode == SEL_Additive )
		{
			double fit = 1;
			for(vectorop::iterator s = m_selectors.begin(), sEnd=m_selectors.end();
				s != sEnd; ++s)
			fit -= 1 - static_cast<selector* >(*s)->indFitness(ind);
			return fit<0?0.:fit;
		}
		else if(m_mode == SEL_Heterogeneity)
			/// fixme
		{
			double fit = 1;
			for(vectorop::iterator s = m_selectors.begin(), sEnd=m_selectors.end();
				s != sEnd; ++s)
			fit *= 1 - static_cast<selector* >(*s)->indFitness(ind);
			return 1-fit;
		}
		// this is the case for SEL_none.
		return 1.0;
	}

	double pySelector::indFitness(individual * ind)
	{
		if( m_len == 0)
		{
			m_len = m_loci.size() * ind->ploidy();
			m_alleles.resize( m_len);
			m_numArray = Allele_Vec_As_NumArray( m_alleles.begin(), m_alleles.end() );
		}

		DBG_FAILIF( static_cast<size_t>(m_len) != ind->ploidy() * m_loci.size(),
			SystemError,
			"Length of m_len is wrong. Have you changed pop type?" );

		UINT pEnd = ind->ploidy();
		for(size_t i=0, iEnd=m_loci.size(), j=0; i < iEnd; ++i)
			for(UINT p=0; p < pEnd; ++p)
				m_alleles[j++] = ind->allele(m_loci[i], p);

		double resDouble;
		PyCallFunc(m_func, "(O)", m_numArray, resDouble, PyObj_As_Double);
		return resDouble;
	}

	PyObject* sample::samples(population& pop)
	{
		DBG_FAILIF( m_name.empty() && m_nameExpr.empty(), ValueError,
			"No sample or sample is not saved. (with empty name and nameExpr)");

		// save sample to local namespace?
		string sampleName="", saveAsName="";
		if( ! m_nameExpr.empty() )
		{
			m_nameExpr.setLocalDict(pop.dict());
			sampleName = m_nameExpr.valueAsString();
		}
		else if( ! m_name.empty() )
			sampleName = m_name;

		try
		{
			PyObject * s = pop.getVar(sampleName);
			Py_INCREF(s);

			return s;
		}
		catch(...)								  // if there is no sample
		{
			Py_INCREF(Py_None);
			return(Py_None);
		}
	}

	bool sample::apply(population& pop)
	{
		// get info from pop
		// if fail, do nothing.
		if(! prepareSample(pop) )
			return true;

		for(UINT t=0; t<m_times; ++t)
		{
			// info of individuals should be set
			population& sample = drawsample(pop);

			// save sample to local namespace?
			string sampleName="", saveAsName="";
			if( ! m_nameExpr.empty() )
			{
				m_nameExpr.setLocalDict(pop.dict());
				sampleName = m_nameExpr.valueAsString();
			}
			else if( ! m_name.empty() )
				sampleName = m_name;

			if(sampleName != "")
			{
				sampleName += "[" + toStr(t) + "]";

				PyObject* popObj = pyPopObj(static_cast<void*>(&sample));

				if( popObj == NULL)
					throw SystemError("Could not expose sample pointer.");

				if( pop.hasVar( sampleName))
					cout << "Warning: sample " << sampleName  <<
						" already exists and will be overwritten."  << endl <<
						"You may want to set name parameter to avoid conflict." << endl;

				// store it.
				// I am not quite sure if this is needed
				Py_INCREF(popObj);
				pop.setVar(sampleName, popObj);
			}

			// save it to a file
			if( ! m_saveAsExpr.empty())
			{
				m_saveAsExpr.setLocalDict(pop.dict());
				saveAsName = m_saveAsExpr.valueAsString();
			}
			else if( ! m_saveAs.empty() )
				saveAsName = m_saveAs;

			if( saveAsName != "" )
			{
				if( m_times > 1 )
					saveAsName = saveAsName + toStr(t);

				sample.savePopulation(saveAsName, m_format);
			}

			if( sampleName == "")
				delete &sample;
		}

		return true;
	}

	bool randomSample::prepareSample(population& pop)
	{
		DBG_FAILIF(m_size.size() == 1 && m_size[0] > pop.popSize(), ValueError,
			"sample size > population size. Can not continue.");

		DBG_FAILIF( m_size.empty(), ValueError,
			"sample size can not be zero/empty for random sampling.");

		DBG_FAILIF( m_size.size() > 1 && m_size.size() != pop.numSubPop(),
			ValueError, "Length of size and number of subpops do not match.");

		if( m_size.size() > 1)
		{
			for(UINT sp = 0; sp < pop.numSubPop(); ++sp)
			{
				DBG_FAILIF( m_size[sp] > pop.subPopSize(sp), ValueError,
					"sample size exceed subpopulation size.");
			}
		}
		return true;
	}

	population& randomSample::drawsample(population& pop)
	{
		pop.setIndSubPopIDWithID();

		if( m_size.size() == 1)
		{
			// randomly choice
			vectori pick(pop.popSize());
			for(size_t i=0, iEnd = pop.popSize(); i< iEnd; ++i)
				pick[i] = i;

			// random shuffle the index array.
			random_shuffle(pick.begin(), pick.end());

			// The first m_size individuals will have its own subpopID

			// remove others
			for(size_t i=m_size[0], iEnd=pop.popSize(); i<iEnd; ++i)
				pop.ind( pick[i] ).setSubPopID( -1 );
		}
		else
		{
			vectori pick;
			for(UINT sp = 0; sp < pop.numSubPop(); ++sp)
			{
				ULONG spSize = pop.subPopSize(sp);
				pick.resize(spSize);
				// randomly choice
				for(size_t i=0, iEnd = spSize; i< iEnd; ++i)
					pick[i] = i;

				// random shuffle the index array.
				random_shuffle(pick.begin(), pick.end());

				// remove others
				for(size_t i=m_size[sp]; i<spSize; ++i)
					pop.ind( pick[i], sp ).setSubPopID( -1 );
			}									  // each subpop
		}										  // else
		return pop.newPopByIndID(-1);
	}

	bool caseControlSample::prepareSample(population& pop )
	{
		if( ! m_spSample )						  // sample from the whole population.
		{
			DBG_FAILIF( m_numCases.size() > 1 || m_numControls.size() > 1,
				ValueError, "Cases, controls need to be a number if sample from the whole population.");

			// first get number of affected.
			int numAffected = count_if(pop.indBegin(), pop.indEnd(),
				isAffected<individual>());
			int numUnaffected = pop.popSize() - numAffected;

			if( m_numCases.size() == 1 && m_numCases[0] > numAffected)
				cout << "Warning: Not enough affected individuals to be sampled: " <<
					"expected " << m_numCases[0] << " available: " << numAffected;

			if( m_numControls.size() == 1 && m_numControls[0] > numUnaffected)
				cout << "Warning: Not enough unaffected individuals to be sampled." <<
					"expected " << m_numControls[0] << " available: " << numUnaffected;

			// save individual indices
			m_caseIdx.resize(1);
			m_controlIdx.resize(1);
			m_caseIdx[0].resize(numAffected);
			m_controlIdx[0].resize(numUnaffected);
			for(size_t i=0, j=0, k=0, iEnd = pop.popSize(); i < iEnd; ++i)
			{
				if( pop.ind(i).affected() )
					m_caseIdx[0][j++] = i;
				else
					m_controlIdx[0][k++] = i;
			}
		}
		else
		{
			UINT numSP = pop.numSubPop();
			vectori numAffected(numSP);
			vectori numUnaffected(numSP);
			m_caseIdx.resize(numSP);
			m_controlIdx.resize(numSP);

			if( m_numCases.size() != numSP || m_numControls.size() != numSP)
				throw ValueError("Size of cases/controls does not match number of subpopulations.");

			for(UINT sp=0; sp < numSP; ++sp)
			{
				// first get number of affected.
				numAffected[sp] = count_if(pop.indBegin(sp), pop.indEnd(sp),
					isAffected<individual>());
				numUnaffected[sp] = pop.subPopSize(sp) - numAffected[sp];

				if( m_numCases[sp] > numAffected[sp])
					cout << "Warning: Not enough affected individuals to be sampled: " <<
						"expected " << m_numCases[sp] << " available: " << numAffected[sp];

				if( m_numControls[sp] > numUnaffected[sp])
					cout << "Warning: Not enough unaffected individuals to be sampled." <<
						"expected " << m_numControls[sp] << " available: " << numUnaffected[sp];

				// save indices
				m_caseIdx[sp].resize(numAffected[sp]);
				m_controlIdx[sp].resize(numUnaffected[sp]);
				for(size_t i=0, j=0, k=0, iEnd = pop.subPopSize(sp); i < iEnd; ++i)
				{
					if( pop.ind(i,sp).affected() )
						m_caseIdx[sp][j++] = i;
					else
						m_controlIdx[sp][k++] = i;
				}
			}
		}
		return true;
	}

	population& caseControlSample::drawsample(population& pop)
	{
		if( ! m_spSample)						  // draw sample from the whole population
		{
			DBG_DO(DBG_SELECTOR, cout << "Selecting from the whole population" << endl);
			// now choose m_caseIdx and m_controlIdx
			// random shuffle the index array.
			random_shuffle(m_caseIdx[0].begin(), m_caseIdx[0].end());
			random_shuffle(m_controlIdx[0].begin(), m_controlIdx[0].end());

			int numAffected = m_caseIdx[0].size();
			int numUnaffected = m_controlIdx[0].size();

			// keep first m_size individuals of shuffled indices
			int nCase, nControl;
			if( m_numCases.empty()  || m_numCases[0] == 0 || m_numCases[0] > numAffected )
				nCase = numAffected;
			else
				nCase = m_numCases[0];

			if( m_numControls.empty()  || m_numControls[0] == 0 || m_numControls[0] > numUnaffected)
				nControl = numUnaffected;
			else
				nControl = m_numControls[0];

			DBG_DO(DBG_SELECTOR, cout << "nCase: " << nCase << " nControl: " << nControl << endl);

			// keep track of which how many case/control from each subpop
			int i;
			vectori nCaseInSP(pop.numSubPop()), nControlInSP(pop.numSubPop());
			for( i=0; i < nCase; ++i)
			{
				nCaseInSP[pop.subPopIndPair(m_caseIdx[0][i]).first]++;
				pop.ind( m_caseIdx[0][i] ).setSubPopID( 0 );
			}
			// remove others
			for( i= nCase; i < numAffected; ++i)
				pop.ind( m_caseIdx[0][i] ).setSubPopID( -1 );

			// keep first m_size individuals of shuffled indices
			for( i=0; i < nControl; ++i)
			{
				nControlInSP[pop.subPopIndPair(m_controlIdx[0][i]).first]++;
				pop.ind( m_controlIdx[0][i] ).setSubPopID( 1 );
			}
			// remove others
			for( i= nControl; i < numUnaffected; ++i)
				pop.ind( m_controlIdx[0][i] ).setSubPopID( -1 );

			DBG_DO(DBG_SELECTOR, cout << "Getting sample population" << endl);
			population& sample = pop.newPopByIndID(-1);

			sample.setIntVectorVar("nCases", nCaseInSP);
			sample.setIntVectorVar("nControls", nControlInSP);
			// determine exactly how many cases and controls in the final sample
			return sample;
		}
		else									  // sample from each subpop
		{
			UINT numSP = pop.numSubPop();
			int nCase, nControl;

			for(UINT sp=0; sp < numSP; ++sp)
			{
				// now choose m_caseIdx and m_controlIdx
				// random shuffle the index array.
				random_shuffle(m_caseIdx[sp].begin(), m_caseIdx[sp].end());
				random_shuffle(m_controlIdx[sp].begin(), m_controlIdx[sp].end());

				// keep first m_size individuals of shuffled indices
				nCase = std::min(m_numCases[sp], static_cast<int>(m_caseIdx[sp].size()));
				for(int i=0; i < nCase; ++i)
					pop.ind( m_caseIdx[sp][i],sp ).setSubPopID( 0 );
				// remove others
				for(int i= nCase, iEnd= m_caseIdx[sp].size(); i<iEnd; ++i)
					pop.ind( m_caseIdx[sp][i],sp ).setSubPopID( -1 );

				// keep first m_size individuals of shuffled indices
				nControl = std::min(m_numControls[sp], static_cast<int>(m_controlIdx[sp].size()));
				for(int i=0; i < nControl; ++i)
					pop.ind( m_controlIdx[sp][i],sp ).setSubPopID( 1 );
				// remove others
				for(int i= nControl, iEnd = m_controlIdx[sp].size(); i<iEnd; ++i)
					pop.ind( m_controlIdx[sp][i],sp ).setSubPopID( -1 );
			}
			// newPop .... but ignore ancestral populations
			population& sample = pop.newPopByIndID(-1);
			sample.setIntVectorVar("nCases", m_numCases);
			sample.setIntVectorVar("nControls", m_numControls);
			return sample;
		}
	}

	bool affectedSibpairSample::prepareSample(population& pop)
	{
		m_father_id = pop.infoIdx(infoField(0));
		m_mother_id = pop.infoIdx(infoField(1));

		// get parental index for each subpop
		DBG_FAILIF( m_size.size() > 1 && m_size.size() != pop.numSubPop(),
			ValueError,
			"Length of array size and number of subpopulations do not match.");

		// sibpairs in each subpopulations
		m_sibpairs.resize(pop.numSubPop());

		size_t nSibs = 0;

		// each subpopulation ...
		for( UINT sp=0; sp < pop.numSubPop(); ++sp)
		{
			// find the sibpairs
			vector< std::pair<ULONG, ULONG> >& sibpairs = m_sibpairs[sp];
			// important!
			sibpairs.clear();

			ULONG spBegin = pop.subPopBegin(sp);
			ULONG spSize = pop.subPopSize(sp);

			// get parents and affected status info
			std::vector< std::pair<ULONG, ULONG> > tags(spSize);
			std::vector<bool> affected(spSize);
			for(size_t i=0; i < spSize; ++i)
			{
				tags[i].first  = static_cast<ULONG>(pop.ind(i, sp).info(m_father_id));
				tags[i].second = static_cast<ULONG>(pop.ind(i, sp).info(m_mother_id));
				affected[i] = pop.ind(i, sp).affected();
			}

			// find affected sibpairs, difficult...
			std::pair<ULONG, ULONG> tag;
			// find an individual with the same parents
			// *** ONLY LOOK AT THE NEXT INDIVIDUAL
			for(size_t i=0, j=1; i< spSize; ++i, ++j)
			{
				if(j == spSize)
					continue;
				if( affected[i] == m_affectedness && affected[j] == m_affectedness
					&& tags[i] == tags[j])
				{
					tag = tags[i];
					// ignore father = mother which is used when no opposite
					// sex exist when mating.
					if(tag.first == tag.second)
						continue;
					// use absolute index, so spBegin+...
					sibpairs.push_back(
						std::pair<ULONG,ULONG>(spBegin+i, spBegin+j) );
					// mark location i, j as occupied.
					affected[i] = ! m_affectedness;
					affected[j] = ! m_affectedness;
					// excluding indiviudus having one of the same parent.
					for(size_t k = j+1; k < spSize; ++k)
					{
						if( tags[k].first == tag.first || tags[k].second == tag.second)
							affected[k] = ! m_affectedness;
					}
					break;
				}
			}									  // for all i

			DBG_DO(DBG_SELECTOR, cout << "Sibpairs at SP " << sp
				<< " is " << sibpairs.size() << endl);

			// now, we know the number of affected sibpairs in subpop i
			pop.setIntVar( subPopVar_String(sp, "numAffectedSibpairs"),
				sibpairs.size());

			nSibs += sibpairs.size();

			if( m_size.size() > 1 && m_size[sp] > sibpairs.size())
				cout << "Warning: Not enough sibpairs (" << sibpairs.size()
					<< ") to be sampled at subpop " << sp << endl;
		}										  // each subpop

		pop.setIntVar( "numAffectedSibpairs", nSibs );

		// do not do sampling if countOnly
		if(m_countOnly)
			return false;
		else
			return true;
	}

	population& affectedSibpairSample::drawsample(population& pop)
	{
		// parents
		vector< std::pair<ULONG, ULONG> > parents;
		// chosen individuals
		vectori nSibpairSP(pop.numSubPop());

		if(m_size.size() <= 1)					  // draw from the whole population
		{
			DBG_DO(DBG_SELECTOR, cout << "Draw from the whole population." << endl);

			/// sibs when m_size.size() <= 1
			vector< std::pair<ULONG, ULONG> > allSibs;

			for(UINT sp =0; sp < pop.numSubPop(); ++sp)
				allSibs.insert( allSibs.end(), m_sibpairs[sp].begin(), m_sibpairs[sp].end());

			if(!m_size.empty() && m_size[0] > allSibs.size())
				cout << "Warning: Not enough sibpairs to be sampled. Requested "
					<< m_size[0] << ", existing " << allSibs.size() << endl;

			ULONG asSize = allSibs.size();

			vectorlu idx(asSize);
			for(size_t i=0; i< asSize; ++i)
				idx[i] = i;

			// set all info to -1
			for(population::IndIterator ind=pop.indBegin();
				ind != pop.indEnd(); ++ind)
			{
				ind->setSubPopID(-1);
			}

			// sample sibpairs
			random_shuffle(idx.begin(), idx.end());

			UINT N=0;
			if(m_size.empty() || m_size[0] > asSize )
				N = asSize;
			else
				N = m_size[0];

			DBG_DO(DBG_SELECTOR, cout << "Getting " << N << " sibpairs" << endl);

			// now we need the first N in order
			std::sort(idx.begin(), idx.begin()+N);

			// set sibpairs
			for( size_t i=0; i< N; ++i)
			{
				pop.ind( allSibs[ idx[i] ].first ).setSubPopID(i);
				pop.ind( allSibs[ idx[i] ].second ).setSubPopID(i);
				std::pair<ULONG, ULONG> tag;
				tag.first = static_cast<ULONG>(pop.ind( allSibs[ idx[i]].first ).info(m_father_id));
				tag.second= static_cast<ULONG>(pop.ind( allSibs[ idx[i]].first ).info(m_mother_id));
				parents.push_back(
					std::pair<ULONG, ULONG>((ULONG)(tag.first), (ULONG)(tag.second))
					);
				DBG_DO(DBG_SELECTOR, cout << i << ": Off: " << int(allSibs[ idx[i] ].first) << " " <<
					int(allSibs[ idx[i] ].second) << " Par: " << pop.ind(allSibs[ idx[i]].first ).info(m_father_id) <<
					" " << pop.ind(allSibs[ idx[i]].first ).info(m_mother_id) << endl);
			}

			// keep the order and count number of selected ones in each subpop
			// here we say number of individuals chosen from a sp += 2
			for( size_t i=0; i<N; ++i)
				nSibpairSP[ pop.subPopIndPair( allSibs[idx[i]].first).first] += 2;
			//
			DBG_DO(DBG_SELECTOR, cout << "#ind from population " << nSibpairSP << endl);
		}
		else									  // for each subpop
		{
			ULONG sibID = 0;
			for(UINT sp = 0; sp < pop.numSubPop(); ++sp)
			{
				vector< std::pair<ULONG, ULONG> >& sibpairs = m_sibpairs[sp];

				size_t asSize = sibpairs.size();
				vectorlu idx(asSize);
				for(size_t i=0; i< asSize; ++i)
					idx[i] = i;

				// set all indices to -1
				for(population::IndIterator ind=pop.indBegin(sp); ind != pop.indEnd(sp); ++ind)
					ind->setSubPopID(-1);

				UINT N = std::min(asSize, static_cast<size_t>(m_size[sp]));

				// sample sibpairs
				random_shuffle(idx.begin(), idx.end());
				std::sort( idx.begin(), idx.begin() + N);

				// set sibpairs
				for( size_t i=0; i< N; ++i)
				{
					// sibpairs uses absolute indices, so no sp parameter is rquired.
					pop.ind( sibpairs[ idx[i]].first ).setSubPopID(sibID);
					pop.ind( sibpairs[ idx[i]].second ).setSubPopID(sibID);
					std::pair<ULONG, ULONG> tag;
					tag.first = pop.ind(sibpairs[ idx[i]].first ).info(m_father_id);
					tag.second= pop.ind(sibpairs[ idx[i]].first ).info(m_mother_id);
					parents.push_back(
						std::pair<ULONG, ULONG>((ULONG)(tag.first), (ULONG)(tag.second))
						);
					sibID ++;
				}
				nSibpairSP[sp] = N*2;
			}									  // sp
		}
		pop.useAncestralPop(1);

		DBG_DO(DBG_SELECTOR, cout << "Working on parent generations." << endl);

		// set everyone to -1 (remove them)
		for(population::IndIterator ind=pop.indBegin();
			ind != pop.indEnd(); ++ind)
		{
			ind->setSubPopID(-1);
		}										  // else
		for(size_t i=0; i< parents.size(); ++i)
		{
			DBG_ASSERT(pop.ind(parents[i].first).subPopID() == -1, SystemError,
				"Duplicate parents are selected: " + toStr(i) + \
				" This can be caused by uni-sex mating so that one individual can be both " + \
				" father and mother ");
			DBG_ASSERT(pop.ind(parents[i].second).subPopID() == -1, SystemError,
				"Duplicate parents are selected: " + toStr(i) + \
				" This can be caused by uni-sex mating so that one individual can be both " + \
				" father and mother ");
			pop.ind(parents[i].first).setSubPopID(i);
			pop.ind(parents[i].second).setSubPopID(i);
		}
		// this is offspring population with copy of ancestral pop
		// (1 means keep only one ancestral population
		population & newPop = pop.newPopByIndID(1);

		DBG_DO(DBG_SELECTOR, cout << "Offspring selection done." << endl);
		return newPop;
	}

	bool largePedigreeSample::prepareSample(population& pop)
	{
		DBG_FAILIF(pop.ancestralDepth() < 2, ValueError,
			"At least two ancestral populations are needed to draw large pedigrees");

		//
		vectorstr fields(2 + m_maxOffspring);
		vectori offspringIdx(m_maxOffspring);
		fields[0] = "pedindex";
		fields[1] = "spouse";
		for(size_t i = 0; i < m_maxOffspring; ++i)
			fields[i+2] = "offspring" + toStr(i);
		// add info fields
		pop.addInfoFields(fields, -1);
		UINT pedindexIdx = pop.infoIdx("pedindex");
		UINT spouseIdx = pop.infoIdx("spouse");
		UINT fatherIdx = pop.infoIdx("father_idx");
		UINT motherIdx = pop.infoIdx("mother_idx");

		for(size_t i = 0; i < m_maxOffspring; ++i)
			offspringIdx[i] = pop.infoIdx(fields[i+2]);
		//
		DBG_DO(DBG_SELECTOR, cout << "Finding spouse and offspring of all individuals" << endl);
		for(size_t ans = 1; ans <=2; ++ans)
		{
			pop.useAncestralPop(ans - 1);
			vectorf dad = pop.indInfo(fatherIdx, true);
			vectorf mom = pop.indInfo(motherIdx, true);
			//
			pop.useAncestralPop(ans);
			//
			for(size_t idx = 0; idx < dad.size(); ++idx)
			{
				double dadSpouse = pop.ind(static_cast<ULONG>(dad[idx])).info(spouseIdx);
				double momSpouse = pop.ind(static_cast<ULONG>(mom[idx])).info(spouseIdx);
				// new case.
				if(dadSpouse == -1. && momSpouse == -1)
				{
					pop.ind(static_cast<ULONG>(dad[idx])).setInfo(mom[idx], spouseIdx);
					pop.ind(static_cast<ULONG>(mom[idx])).setInfo(dad[idx], spouseIdx);
					pop.ind(static_cast<ULONG>(dad[idx])).setInfo(idx, offspringIdx[0]);
					pop.ind(static_cast<ULONG>(mom[idx])).setInfo(idx, offspringIdx[0]);
				}
				else if(dadSpouse != -1 && momSpouse != -1 &&
					dadSpouse == mom[idx] && momSpouse == dad[idx])
				{
					// which offspring
					for(size_t i=1; i < m_maxOffspring; ++i)
					{
						if(pop.ind(dad[idx]).info(offspringIdx[i]) == -1)
						{
							pop.ind(static_cast<ULONG>(dad[idx])).setInfo(idx, offspringIdx[i]);
							pop.ind(static_cast<ULONG>(mom[idx])).setInfo(idx, offspringIdx[i]);
						}
					}
				}
			}
		}										  // ans = 1, 2
		//
		pop.useAncestralPop(2);
		size_t g3size = pop.popSize();
		size_t pedindex = 1;
		//
		DBG_DO(DBG_SELECTOR, cout << "Finding all three-generation pedigrees" << endl);
		for(size_t idx=0; idx < g3size; ++idx)
		{
			pop.useAncestralPop(2);
			int pedSize = 2;
			int numAffected = 0;
			//
			// already belong to other pedigree
			int grandspouse = static_cast<int>(pop.ind(idx).info(spouseIdx));
			// no spuse? one of grandparents belong to another pedigree?
			if (grandspouse < 0 || pop.ind(idx).info(pedindexIdx) >= 1 ||
				pop.ind(grandspouse).info(pedindexIdx) >= 1)
				continue;
			if(pop.ind(idx).affected())
				numAffected ++;
			if(pop.ind(grandspouse).affected())
				numAffected ++;
			//
			vectorlu parentsTmp;
			for(size_t x = 0; x < m_maxOffspring; ++x)
			{
				int off = static_cast<int>(pop.ind(idx).info(offspringIdx[x]));
				if (off != -1)
					parentsTmp.push_back(off);
			}
			// parents generation
			//
			if(parentsTmp.empty())
				continue;
			// go to parent generation
			pop.useAncestralPop(1);
			vectorlu spouseofparents;
			vectorlu childrenTmp;
			// verify parents
			vectorlu parents;
			for(vectorlu::iterator par = parentsTmp.begin(); par != parentsTmp.end(); ++par)
			{
				if(pop.ind(*par).info(pedindexIdx) == -1.)
				{
					parents.push_back(*par);
					pedSize ++;
					if(pop.ind(*par).affected())
						numAffected ++;
				}
			}
			//
			for(vectorlu::iterator par = parents.begin(); par != parents.end(); ++par)
			{
				int spouse = static_cast<int>(pop.ind(*par).info(spouseIdx));
				// if there is spouse, add it in
				if(spouse >= 0)
				{
					spouseofparents.push_back(spouse);
					pedSize ++;
					if(pop.ind(spouse).affected())
						numAffected ++;
					// there are children only when there is spouse
					for(size_t x = 0; x < m_maxOffspring; ++x)
					{
						int off = static_cast<int>(pop.ind(*par).info(offspringIdx[x]));
						if (off != -1)
							childrenTmp.push_back(off);
					}
				}
			}
			if(childrenTmp.empty())
				continue;
			pop.useAncestralPop(0);
			// verify children
			vectorlu children;
			for(vectorlu::iterator child = childrenTmp.begin(); child != childrenTmp.end(); ++child)
			{
				if(pop.ind(*child).info(pedindexIdx) == -1.)
				{
					children.push_back(*child);
					pedSize ++;
					if(pop.ind(*child).affected())
						numAffected ++;
				}
			}
			if(children.empty())
				continue;
			if(pedSize < m_minPedSize || (m_minAffected > 0 && numAffected < m_minAffected))
				continue;
			// everything seems to be fine, add them to the family
			pop.useAncestralPop(2);
			// grandparents
			pop.ind(idx).setInfo(pedindex, pedindexIdx);
			pop.ind(grandspouse).setInfo(pedindex, pedindexIdx);
			// parents and their spouse
			pop.useAncestralPop(1);
			for(vectorlu::iterator it=parents.begin(); it != parents.end(); ++it)
				pop.ind(*it).setInfo(pedindex, pedindexIdx);
			// spouse of parents
			for(vectorlu::iterator it=spouseofparents.begin(); it != spouseofparents.end(); ++it)
				pop.ind(*it).setInfo(pedindex, pedindexIdx);
			// now children
			pop.useAncestralPop(0);
			for(vectorlu::iterator it=children.begin(); it != children.end(); ++it)
				pop.ind(*it).setInfo(pedindex, pedindexIdx);
			// is this family qualified?
			m_validPedigrees.push_back(pedindex);
			pedindex += 1;
		}
		pop.useAncestralPop(0);
		// print the pedigrees
		if(m_validPedigrees.size() < m_size)
			cout << "Only " << m_validPedigrees.size() << "valid pedigrees are found. " << endl;

		// do not do sampling if countOnly
		if(m_countOnly)
			return false;
		else
			return true;
	}

	population& largePedigreeSample::drawsample(population& pop)
	{
		// sample sibpairs
		DBG_DO(DBG_SELECTOR, cout << "Generating sample" << endl);
		size_t N = m_size;
		if(m_validPedigrees.size() > m_size)
			random_shuffle(m_validPedigrees.begin(), m_validPedigrees.end());
		else
			N = m_validPedigrees.size();

		DBG_DO(DBG_SELECTOR, cout << "Getting " << N << " pedigrees" << endl);

		// now we need the first N in order
		std::sort(m_validPedigrees.begin(), m_validPedigrees.begin()+N);

		// set sibpairs
		size_t spouseIdx = pop.infoIdx("spouse");
		size_t pedindexIdx = pop.infoIdx("pedindex");
		vectorf grandIdx = pop.indInfo(pedindexIdx, true);
		for(vectori::iterator ped=m_validPedigrees.begin(); ped != m_validPedigrees.begin()+N; ++ped)
		{
			// find these offspring....
			pop.useAncestralPop(2);
			for(size_t i = 0; i<pop.popSize(); ++i)
				pop.ind(i).setSubPopID(-1);
			// find grandparent
			vectorf::iterator tmp = find(grandIdx.begin(), grandIdx.end(), *ped);
			DBG_FAILIF(tmp==grandIdx.end(), ValueError, "Can not find pedigree");
			//
			vectori offspringIdx(m_maxOffspring);
			for(size_t i=0; i<m_maxOffspring; ++i)
				offspringIdx[i] = pop.infoIdx("offspring"+toStr(i));
			size_t grandpar1 = tmp - grandIdx.begin();
			size_t grandpar2 = static_cast<size_t>(pop.ind(grandpar1).info(spouseIdx));
			pop.ind(grandpar1).setSubPopID(*ped);
			pop.ind(grandpar2).setSubPopID(*ped);
			// find parents
			vectorlu parents;
			for(size_t x = 0; x < m_maxOffspring; ++x)
			{
				size_t off = static_cast<size_t>(pop.ind(grandpar1).info(offspringIdx[x]));
				if(off != -1)
					parents.push_back(off);
			}
			//
			pop.useAncestralPop(1);
			vectorlu children;
			for(vectorlu::iterator it=parents.begin(); it != parents.end(); ++it)
			{
				if(pop.ind(*it).info(pedindexIdx) == *ped)
				{
					pop.ind(*it).setSubPopID(*ped);
					size_t spouse =  static_cast<size_t>(pop.ind(*it).info(spouseIdx));
					// if there is spouse, add it in
					if(spouse >= 0)
					{
						pop.ind(spouse).setSubPopID(*ped);
						for(size_t x = 0; x < m_maxOffspring; ++x)
						{
							size_t off = static_cast<size_t>(pop.ind(grandpar1).info(offspringIdx[x]));
							if(off != -1)
								children.push_back(off);
						}
					}
				}
			}
			// go to children
			pop.useAncestralPop(0);
			for(vectorlu::iterator it=children.begin(); it != children.end(); ++it)
			{
				if(pop.ind(*it).info(pedindexIdx) == *ped)
					pop.ind(*it).setSubPopID(*ped);
			}
		}
		population & newPop = pop.newPopByIndID(2);

		DBG_DO(DBG_SELECTOR, cout << "Offspring selection done." << endl);
		return newPop;
	}

}
