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

#include "sampler.h"

namespace simuPOP
{
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

	void sample::saveIndIndex(population& pop, const string& indexField)
	{
		UINT fieldIdx = pop.infoIdx(indexField);
		for(size_t ans = 0; ans <= pop.ancestralDepth(); ++ans)
		{
			pop.useAncestralPop(ans);
			// need to save old index
			for(size_t idx = 0; idx < pop.popSize(); ++idx)
				pop.ind(idx).setInfo(idx, fieldIdx);
		}
		pop.useAncestralPop(0);
	}

	void sample::resetParentalIndex(population& pop, const string& fatherField,
		const string& motherField, const string& indexField)
	{
		UINT fatherIdx = pop.infoIdx(fatherField);
		UINT motherIdx = pop.infoIdx(motherField);
		UINT indexIdx = pop.infoIdx(indexField);

		// for the top generation, no parents
		pop.useAncestralPop(pop.ancestralDepth());
		for(population::IndIterator it=pop.indBegin(); it != pop.indEnd(); ++it)
		{
			// -1 means no parents.
			// we do not use 0 since 0 is valid
			it->setInfo(-1, fatherIdx);
			it->setInfo(-1, motherIdx);
		}
		// for other generations
		for(size_t ans = 0; ans < pop.ancestralDepth(); ++ans)
		{
			// parents... get their old index
			pop.useAncestralPop(ans+1);
			vectorf oldindex = pop.indInfo(indexIdx, true);
			// children
			pop.useAncestralPop(ans);
			for(population::IndIterator it=pop.indBegin(); it != pop.indEnd(); ++it)
			{
				// what is the idx of my parents now?
				vectorf::iterator tmp = find(oldindex.begin(), oldindex.end(), it->info(fatherIdx));
				// no parents
				if(tmp == oldindex.end())
					it->setInfo(-1, fatherIdx);
				else
					it->setInfo(tmp - oldindex.begin(), fatherIdx);
				tmp = find(oldindex.begin(), oldindex.end(), it->info(motherIdx));
				if(tmp == oldindex.end())
					it->setInfo(-1, motherIdx);
				else
					it->setInfo(tmp-oldindex.begin(), motherIdx);
			}
		}
		pop.useAncestralPop(0);
	}

	void sample::findOffspringAndSpouse(population& pop, unsigned ancestralDepth,
		unsigned maxOffspring,
		const string& fatherField, const string& motherField,
		const string& spouseField, const string& offspringField)
	{
		vectori offspringIdx(maxOffspring);
		for(size_t i = 0; i < maxOffspring; ++i)
			offspringIdx[i] = pop.infoIdx(offspringField + toStr(i));
		UINT spouseIdx = pop.infoIdx(spouseField);
		UINT fatherIdx = pop.infoIdx(fatherField);
		UINT motherIdx = pop.infoIdx(motherField);

		DBG_DO(DBG_SELECTOR, cout << "Finding spouse and offspring of all individuals" << endl);
		DBG_FAILIF(ancestralDepth > pop.ancestralDepth(), ValueError,
			"The population does not have enough ancestral generations for this operation");

		for(unsigned ans = 1; ans <= ancestralDepth; ++ans)
		{
			pop.useAncestralPop(ans - 1);
			vectorf dad = pop.indInfo(fatherIdx, true);
			vectorf mom = pop.indInfo(motherIdx, true);
			//
			pop.useAncestralPop(ans);
			//
			for(size_t idx = 0; idx < dad.size(); ++idx)
			{
				InfoType dadSpouse = pop.ind(static_cast<ULONG>(dad[idx])).info(spouseIdx);
				InfoType momSpouse = pop.ind(static_cast<ULONG>(mom[idx])).info(spouseIdx);
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
					for(size_t i=1; i < maxOffspring; ++i)
					{
						if(pop.ind(static_cast<ULONG>(dad[idx])).info(offspringIdx[i]) == -1.)
						{
							pop.ind(static_cast<ULONG>(dad[idx])).setInfo(idx, offspringIdx[i]);
							pop.ind(static_cast<ULONG>(mom[idx])).setInfo(idx, offspringIdx[i]);
							break;
						}
					}
				}
			}									  // idx
		}										  // ancestal generations
		pop.useAncestralPop(0);
	}

	void sample::resetSubPopID(population& pop)
	{
		int oldGen = pop.ancestralGen();
		for(size_t anc = 0; anc <= pop.ancestralDepth(); ++anc)
		{
			pop.useAncestralPop(anc);
			for(population::IndIterator it=pop.indBegin(); it != pop.indEnd(); ++it)
				it->setSubPopID(-1);
		}
		pop.useAncestralPop(oldGen);
	}

	bool randomSample::prepareSample(population& pop)
	{
		DBG_FAILIF(m_size.size() == 1 && m_size[0] > pop.popSize(), ValueError,
			"sample size > population size. Can not continue.");

		DBG_FAILIF( m_size.empty(), ValueError,
			"sample size can not be zero/empty for random sampling.");

		DBG_FAILIF( m_size.size() > 1 && m_size.size() != pop.numSubPop(),
			ValueError, "Length of size and number of subpops do not match.");

		pop.addInfoField("oldindex", -1);
		saveIndIndex(pop, "oldindex");
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
		// record the index in the old population, for samples.
		pop.addInfoField("oldindex", -1);
		saveIndIndex(pop, "oldindex");
		if(!m_spSample)							  // sample from the whole population.
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
			population& sample = pop.newPopByIndID(0);

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
			population& sample = pop.newPopByIndID(0);
			sample.setIntVectorVar("nCases", m_numCases);
			sample.setIntVectorVar("nControls", m_numControls);
			return sample;
		}
	}

	bool affectedSibpairSample::prepareSample(population& pop)
	{
		// get parental index for each subpop
		DBG_FAILIF( m_size.size() > 1 && m_size.size() != pop.numSubPop(),
			ValueError,
			"Length of array size and number of subpopulations do not match.");

		m_validSibs.clear();

		m_father_id = pop.infoIdx(infoField(0));
		m_mother_id = pop.infoIdx(infoField(1));
		vectorstr fields(5);
		fields[0] = "oldindex";
		fields[1] = "pedindex";
		fields[2] = "offspring0";
		fields[3] = "offspring1";
		fields[4] = "spouse";
		pop.addInfoFields(fields, -1);
		saveIndIndex(pop, "oldindex");
		UINT pedindexIdx = pop.infoIdx("pedindex");
		UINT off0Idx = pop.infoIdx("offspring0");
		UINT off1Idx = pop.infoIdx("offspring1");
		UINT spouseIdx = pop.infoIdx("spouse");

		// 1: one ancestralDepth
		// 2: two offsprings
		findOffspringAndSpouse(pop, 1, 2, "father_idx", "mother_idx",
			"spouse", "offspring");				  // ans = 1, 2
		//
		// find sibpairs from the parental generation.
		pop.useAncestralPop(1);
		vectorlu off;
		int pedIdx = 0;
		// valid sibpairs for each subpopulation
		for( UINT sp = 0; sp < pop.numSubPop(); ++sp)
		{
			for(population::IndIterator it = pop.indBegin(sp); it != pop.indEnd(sp); ++it)
			{
				// individual already belongs to another family
				if(it->info(pedindexIdx) != -1.)
					continue;
				double spouse = it->info(spouseIdx);
				// has spouse, spouse does not belong to anther ped, and two kids
				if(spouse != -1. && pop.ind(static_cast<UINT>(spouse)).info(pedindexIdx) == -1.
					&& it->info(off1Idx) != -1.)
				{
					it->setInfo(pedIdx, pedindexIdx);
					pop.ind(static_cast<UINT>(spouse)).setInfo(pedIdx, pedindexIdx);
					off.push_back(static_cast<ULONG>(it->info(off0Idx)));
					off.push_back(static_cast<ULONG>(it->info(off1Idx)));
					pedIdx++;
				}
			}
			DBG_DO(DBG_SELECTOR, cout << "Number of sibpairs in subpop " << sp << " is "
				<< m_validSibs[sp].size() << endl);
		}										  // each subpop
		pop.useAncestralPop(0);
		m_validSibs.resize(pop.numSubPop());
		for( UINT sp = 0; sp < pop.numSubPop(); ++sp)
			m_validSibs[sp].clear();
		UINT nSibs = 0;
		for(size_t i = 0; i < off.size()/2; ++i)
		{
			pop.ind(off[2*i]).setInfo(i, pedindexIdx);
			pop.ind(off[2*i+1]).setInfo(i, pedindexIdx);
			if( m_affectedness == pop.ind(off[2*i]).affected() 
				&& m_affectedness == pop.ind(off[2*i+1]).affected())
			{
				m_validSibs[pop.subPopIndPair(off[2*i]).first].push_back(i);
				nSibs ++;
			}
		}
		pop.setIntVar("numAffectedSibpairs", nSibs);
		// do not do sampling if countOnly
		if(m_countOnly)
			return false;
		else
			return true;
	}

	population& affectedSibpairSample::drawsample(population& pop)
	{
		// mark to remove everyone
		resetSubPopID(pop);
		vectorlu acceptedSibs;

		if(m_size.size() <= 1)					  // draw from the whole population
		{
			/// collect all families
			vector<UINT> allSibs;

			for(UINT sp =0; sp < pop.numSubPop(); ++sp)
				allSibs.insert(allSibs.end(), m_validSibs[sp].begin(), m_validSibs[sp].end());

			if(!m_size.empty() && m_size[0] > allSibs.size())
				cout << "Warning: Not enough sibpairs to be sampled. Requested "
					<< m_size[0] << ", existing " << allSibs.size() << endl;

			UINT N = 0;
			if(m_size.empty() || m_size[0] > allSibs.size())
				N = allSibs.size();
			else
			{
				N = m_size[0];
				// random shuffle.
				random_shuffle(allSibs.begin(), allSibs.end());
			}
			acceptedSibs.insert(acceptedSibs.end(), allSibs.begin(), allSibs.begin() + N);
		}
		else									  // for each subpop
		{
			for(UINT sp = 0; sp < pop.numSubPop(); ++sp)
			{
				vectorlu & sibpairs = m_validSibs[sp];

				UINT N = sibpairs.size();
				if(N > m_size[sp])
				{
					N = m_size[sp];
					// sample sibpairs
					random_shuffle(sibpairs.begin(), sibpairs.end());
				}
				acceptedSibs.insert(acceptedSibs.end(), sibpairs.begin(), sibpairs.begin() + N);
			}									  // sp
		}
		// now, we have acepted sibs, set subpopid in preparation for a new
		// population
		pop.useAncestralPop(1);

		vectorlu offspring;
		int pedIdx = 0;
		vectorlu off;
		UINT spouseIdx = pop.infoIdx("spouse");
		UINT pedindexIdx = pop.infoIdx("pedindex");
		UINT off0Idx = pop.infoIdx("offspring0");
		UINT off1Idx = pop.infoIdx("offspring1");

		for(size_t i = 0; i < pop.popSize(); ++i)
		{
			individual& ind = pop.ind(i);
			double infoPedIdx = ind.info(pedindexIdx);
			if(infoPedIdx == -1.)
				continue;
			double spouse = ind.info(spouseIdx);
			// only look forward
			if(spouse < i)
				continue;
			// if this family is selected
			if(find(acceptedSibs.begin(), acceptedSibs.end(), static_cast<ULONG>(infoPedIdx)) != acceptedSibs.end())
			{
				ind.setSubPopID(pedIdx);
				pop.ind(static_cast<UINT>(spouse)).setSubPopID(pedIdx);
				off.push_back(static_cast<ULONG>(ind.info(off0Idx)));
				off.push_back(static_cast<ULONG>(ind.info(off1Idx)));
				pedIdx ++;
			}
		}
		pop.useAncestralPop(0);
		for(size_t i = 0; i < off.size()/2; ++i)
		{
			pop.ind(off[2*i]).setSubPopID(i);
			pop.ind(off[2*i+1]).setSubPopID(i);
		}

		// this is offspring population with copy of ancestral pop
		// (1 means keep only one ancestral population
		population & newPop = pop.newPopByIndID(1);
		//
		resetParentalIndex(newPop, "father_idx", "mother_idx", "oldindex");

		DBG_DO(DBG_SELECTOR, cout << "Offspring selection done." << endl);
		return newPop;
	}

	bool largePedigreeSample::prepareSample(population& pop)
	{
		DBG_FAILIF(pop.ancestralDepth() < 2, ValueError,
			"At least two ancestral populations are needed to draw large pedigrees");
		//
		m_validPedigrees.clear();
		UINT nPed = 0;
		//
		vectorstr fields(3 + m_maxOffspring);
		vectori offspringIdx(m_maxOffspring);
		for(size_t i = 0; i < m_maxOffspring; ++i)
			fields[i] = "offspring" + toStr(i);
		fields[m_maxOffspring] = "pedindex";
		fields[m_maxOffspring+1] = "spouse";
		fields[m_maxOffspring+2] = "oldindex";
		// add info fields
		pop.addInfoFields(fields, -1);
		UINT pedindexIdx = pop.infoIdx("pedindex");
		UINT spouseIdx = pop.infoIdx("spouse");
		for(size_t i = 0; i < m_maxOffspring; ++i)
			offspringIdx[i] = pop.infoIdx(fields[i]);
		// save old index
		saveIndIndex(pop, "oldindex");
		//
		// 2 means find till grandfather.
		findOffspringAndSpouse(pop, 2, m_maxOffspring, "father_idx", "mother_idx",
			"spouse", "offspring");				  // ans = 1, 2
		//
		pop.useAncestralPop(2);
		m_validPedigrees.resize(pop.numSubPop());
		size_t pedindex = 0;
		DBG_DO(DBG_SELECTOR, cout << "Finding all three-generation pedigrees" << endl);
		for(UINT sp = 0; sp < pop.numSubPop(); ++sp)
		{
			m_validPedigrees[sp].clear();
			size_t g3start = pop.subPopBegin(sp);
			size_t g3end = pop.subPopEnd(sp);
			//
			for(size_t idx = g3start; idx < g3end; ++idx)
			{
				pop.useAncestralPop(2);
				unsigned pedSize = 2;
				unsigned numAffected = 0;
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
					InfoType spouse = pop.ind(*par).info(spouseIdx);
					// if there is spouse, add it in
					if(spouse >= 0. && pop.ind(*par).info(pedindexIdx) == -1.)
					{
						spouseofparents.push_back(static_cast<ULONG>(spouse));
						pedSize ++;
						if(pop.ind(static_cast<ULONG>(spouse)).affected())
							numAffected ++;
						// there are children only when there is spouse
						for(size_t x = 0; x < m_maxOffspring; ++x)
						{
							InfoType off = pop.ind(*par).info(offspringIdx[x]);
							if (off != -1.)
								childrenTmp.push_back(static_cast<ULONG>(off));
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
					// unoccupied.
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
				m_validPedigrees[sp].push_back(boost::tie(pedindex, pedSize));
				pedindex += 1;
			}
			nPed += m_validPedigrees[sp].size();
		}
		pop.useAncestralPop(0);
		pop.setIntVar("numPedigrees", nPed);

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
		resetSubPopID(pop);
		pedArray acceptedPeds;

		if(m_size.size() <= 1)					  // draw from the whole population
		{
			/// collect all families
			pedArray allPeds;

			for(UINT sp =0; sp < pop.numSubPop(); ++sp)
				allPeds.insert(allPeds.end(), m_validPedigrees[sp].begin(), m_validPedigrees[sp].end());

			if(!m_size.empty() && m_size[0] > allPeds.size())
				cout << "Warning: Not enough sibpairs to be sampled. Requested "
					<< m_size[0] << ", existing " << allPeds.size() << endl;

			random_shuffle(allPeds.begin(), allPeds.end());
			UINT N = 0;
			// consider total size.
			if(m_size.empty() || m_minTotalSize > 0)
			{
				size_t totalSize = 0;
				for(pedArray::iterator it = allPeds.begin(); it != allPeds.end(); ++it)
				{
					totalSize += boost::get<1>(*it);
					N++;
					if(totalSize > m_minTotalSize)
						break;
				}
				if(totalSize < m_minTotalSize)
					cout << "Warning: can not reach min total size" << endl;
			}
			else
				N = m_size[0];
			acceptedPeds.insert(acceptedPeds.end(), allPeds.begin(), allPeds.begin() + N);
		}
		else									  // for each subpop
		{
			for(UINT sp = 0; sp < pop.numSubPop(); ++sp)
			{
				pedArray & peds = m_validPedigrees[sp];

				UINT N = peds.size();
				if(N > m_size[sp])
				{
					N = m_size[sp];
					// sample sibpairs
					random_shuffle(peds.begin(), peds.end());
				}
				acceptedPeds.insert(acceptedPeds.end(), peds.begin(), peds.begin() + N);
			}									  // sp
		}

		DBG_DO(DBG_SELECTOR, cout << "Sampling " << acceptedPeds.size() << " pedigrees " << endl);

		//
		UINT spouseIdx = pop.infoIdx("spouse");
		UINT pedindexIdx = pop.infoIdx("pedindex");
		vectori offspringIdx(m_maxOffspring);
		for(size_t i=0; i<m_maxOffspring; ++i)
			offspringIdx[i] = pop.infoIdx("offspring" + toStr(i));
		pop.useAncestralPop(2);
		vectorf grandIdx = pop.indInfo(pedindexIdx, true);
		int newPedID = 0;
		for(pedArray::iterator ped = acceptedPeds.begin(); ped != acceptedPeds.end(); ++ped, ++newPedID)
		{
			double pedID = boost::get<0>(*ped);
			DBG_DO(DBG_SELECTOR, cout << "Getting pedigree " << pedID << " of size " << boost::get<1>(*ped) << endl);
			// real pedsize, for verification purpose
			int ps = 0;
			// find these offspring....
			pop.useAncestralPop(2);
			// find grandparent
			vectorf::iterator tmp = find(grandIdx.begin(), grandIdx.end(), pedID);
			DBG_FAILIF(tmp==grandIdx.end(), ValueError, "Can not find pedigree");
			//
			size_t grandpar1 = tmp - grandIdx.begin();
			DBG_FAILIF(pop.ind(grandpar1).info(spouseIdx) == -1., SystemError,
				"Grand parent's spouse is invalid");
			size_t grandpar2 = static_cast<size_t>(pop.ind(grandpar1).info(spouseIdx));
			pop.ind(grandpar1).setSubPopID(newPedID);
			pop.ind(grandpar2).setSubPopID(newPedID);
			ps += 2;
			// find parents
			vectorlu parents;
			for(size_t x = 0; x < m_maxOffspring; ++x)
			{
				InfoType off = pop.ind(grandpar1).info(offspringIdx[x]);
				if(off != -1.)
					parents.push_back(static_cast<ULONG>(off));
			}
			//
			pop.useAncestralPop(1);
			vectorlu children;
			for(vectorlu::iterator it=parents.begin(); it != parents.end(); ++it)
			{
				if(pop.ind(*it).info(pedindexIdx) == pedID)
				{
					pop.ind(*it).setSubPopID(newPedID);
					ps ++;
					InfoType spouse = pop.ind(*it).info(spouseIdx);
					// if there is spouse, add it in
					if(spouse != -1)
					{
						pop.ind(static_cast<ULONG>(spouse)).setSubPopID(newPedID);
						ps ++;
						DBG_ASSERT(pop.ind(static_cast<ULONG>(spouse)).info(pedindexIdx) == pedID,
							ValueError, "Spouse does not belong to this pedigree, something wrong.");
						for(size_t x = 0; x < m_maxOffspring; ++x)
						{
							InfoType off = pop.ind(*it).info(offspringIdx[x]);
							if(off != -1)
								children.push_back(static_cast<ULONG>(off));
						}
					}
				}
			}
			// go to children
			pop.useAncestralPop(0);
			for(vectorlu::iterator it=children.begin(); it != children.end(); ++it)
			{
				if(pop.ind(*it).info(pedindexIdx) == pedID)
				{
					pop.ind(*it).setSubPopID(newPedID);
					ps ++;
				}
			}
			DBG_FAILIF(ps != boost::get<1>(*ped), SystemError,
				"Pedigree sizes do not match, estimated " + toStr(boost::get<1>(*ped)) + " real: " + toStr(ps));
		}
		// just to make sure
		pop.useAncestralPop(0);
		// saving samples in a new population, 
		population & newPop = pop.newPopByIndID(2);
		//
		resetParentalIndex(newPop, "father_idx", "mother_idx", "oldindex");
		return newPop;
	}
	
	bool nuclearFamilySample::prepareSample(population& pop)
	{
		DBG_FAILIF(pop.ancestralDepth() < 1, ValueError,
			"At least one ancestral populations are needed to draw large pedigrees");
		//
		m_validPedigrees.clear();
		UINT nPed = 0;
		//
		vectorstr fields(3 + m_maxOffspring);
		vectori offspringIdx(m_maxOffspring);
		for(size_t i = 0; i < m_maxOffspring; ++i)
			fields[i] = "offspring" + toStr(i);
		fields[m_maxOffspring] = "pedindex";
		fields[m_maxOffspring+1] = "spouse";
		fields[m_maxOffspring+2] = "oldindex";
		// add info fields
		pop.addInfoFields(fields, -1);
		UINT pedindexIdx = pop.infoIdx("pedindex");
		UINT spouseIdx = pop.infoIdx("spouse");
		for(size_t i = 0; i < m_maxOffspring; ++i)
			offspringIdx[i] = pop.infoIdx(fields[i]);
		// save old index
		saveIndIndex(pop, "oldindex");
		//
		// 1 means find till parents
		findOffspringAndSpouse(pop, 1, m_maxOffspring, "father_idx", "mother_idx",
			"spouse", "offspring");
		// offspring index
		m_validPedigrees.resize(pop.numSubPop());
		size_t pedIdx = 0;
		DBG_DO(DBG_SELECTOR, cout << "Finding all two-generation pedigrees" << endl);
		for(UINT sp = 0; sp < pop.numSubPop(); ++sp)
		{
			pop.useAncestralPop(1);
			vectorf off;
			for(population::IndIterator it = pop.indBegin(sp); it != pop.indEnd(sp); ++it)
			{
				// individual already belongs to another family
				if(it->info(pedindexIdx) != -1.)
					continue;
				double spouse = it->info(spouseIdx);
				// has spouse, spouse does not belong to anther ped
				if(spouse != -1. && pop.ind(static_cast<UINT>(spouse)).info(pedindexIdx) == -1.)
				{
					it->setInfo(pedIdx, pedindexIdx);
					pop.ind(static_cast<UINT>(spouse)).setInfo(pedIdx, pedindexIdx);
					// many of the offspring may be -1.
					for(UINT oi=0; oi < m_maxOffspring; ++oi)
						off.push_back(static_cast<ULONG>(it->info(offspringIdx[oi])));
					pedIdx++;
				}
			}
			// go to offspring generation and verify
			pop.useAncestralPop(0);
			m_validPedigrees[sp].clear();
			for(UINT i=0; i< off.size()/m_maxOffspring; ++i)
			{
				UINT pedSize = 2;
				UINT pedAffected = 0;
				for(UINT oi=0; oi < m_maxOffspring; ++oi)
				{
					// valid offspring
					if (off[i*m_maxOffspring + oi] != -1)
					{
						pop.ind(static_cast<UINT>(off[i*m_maxOffspring + oi])).setInfo(i, pedindexIdx);
						pedSize ++;
						if (pop.ind(static_cast<UINT>(off[i*m_maxOffspring + oi])).affected())
							pedAffected ++;
					}
				}
				if (pedSize > m_minPedSize && pedAffected > m_minAffected)
					m_validPedigrees[sp].push_back(boost::tie(i, pedSize));
			}
			nPed += m_validPedigrees[sp].size();
			DBG_DO(DBG_SELECTOR, cout << "Number of valid pedigrees in subpop " + toStr(sp) + " is "
				+ toStr(m_validPedigrees[sp].size()) << endl);
		}
		pop.useAncestralPop(0);
		pop.setIntVar("numPedigrees", nPed);

		// do not do sampling if countOnly
		if(m_countOnly)
			return false;
		else
			return true;
	}

	population& nuclearFamilySample::drawsample(population& pop)
	{
		DBG_DO(DBG_SELECTOR, cout << "Generating nuclear family" << endl);
		resetSubPopID(pop);
		pedArray acceptedPeds;

		if(m_size.size() <= 1)					  // draw from the whole population
		{
			/// collect all families
			pedArray allPeds;

			for(UINT sp =0; sp < pop.numSubPop(); ++sp)
				allPeds.insert(allPeds.end(), m_validPedigrees[sp].begin(), m_validPedigrees[sp].end());

			if(!m_size.empty() && m_size[0] > allPeds.size())
				cout << "Warning: Not enough nuclear families to be sampled. Requested "
					<< m_size[0] << ", existing " << allPeds.size() << endl;

			random_shuffle(allPeds.begin(), allPeds.end());
			UINT N = 0;
			// consider total size.
			if(m_size.empty() || m_minTotalSize > 0)
			{
				size_t totalSize = 0;
				for(pedArray::iterator it = allPeds.begin(); it != allPeds.end(); ++it)
				{
					totalSize += boost::get<1>(*it);
					N++;
					if(totalSize > m_minTotalSize)
						break;
				}
				if(totalSize < m_minTotalSize)
					cout << "Warning: can not reach min total size" << endl;
			}
			else
				N = m_size[0];
			acceptedPeds.insert(acceptedPeds.end(), allPeds.begin(), allPeds.begin() + N);
		}
		else									  // for each subpop
		{
			for(UINT sp = 0; sp < pop.numSubPop(); ++sp)
			{
				pedArray & peds = m_validPedigrees[sp];

				UINT N = peds.size();
				if(N > m_size[sp])
				{
					N = m_size[sp];
					random_shuffle(peds.begin(), peds.end());
				}
				acceptedPeds.insert(acceptedPeds.end(), peds.begin(), peds.begin() + N);
			}									  // sp
		}

		DBG_DO(DBG_SELECTOR, cout << "Sampling " << acceptedPeds.size() << " pedigrees " << endl);

		//
		UINT spouseIdx = pop.infoIdx("spouse");
		UINT pedindexIdx = pop.infoIdx("pedindex");
		vectori offspringIdx(m_maxOffspring);
		for(size_t i=0; i<m_maxOffspring; ++i)
			offspringIdx[i] = pop.infoIdx("offspring" + toStr(i));
		pop.useAncestralPop(1);
		vectorf off;
		int pedIdx = 0;
		
		for(size_t i=0; i< pop.popSize(); ++i)
		{
			individual& ind = pop.ind(i);
			double infoPedIdx = ind.info(pedindexIdx);
			if (infoPedIdx == -1.)
				continue;
			double spouse = ind.info(spouseIdx);
			// only look forward
			if (static_cast<size_t>(spouse) < i)
				continue;
			// if this family is selected
			for(pedArray::iterator it = acceptedPeds.begin(); it != acceptedPeds.end(); ++it)
			{
				if ( boost::get<0>(*it) == static_cast<ULONG>(infoPedIdx))
				{
					ind.setSubPopID(pedIdx);
					pop.ind(static_cast<UINT>(spouse)).setSubPopID(pedIdx);
					for(UINT oi=0; oi < m_maxOffspring; ++oi)
						off.push_back(ind.info(offspringIdx[oi]));
					pedIdx ++;
					break;
				}
			}
		}
		pop.useAncestralPop(0);
		for(size_t i=0; i<off.size(); ++i)
		{
			if(off[i] != -1.)
				pop.ind(static_cast<ULONG>(off[i])).setSubPopID(i/m_maxOffspring);
		}
		population & newPop = pop.newPopByIndID(1);
		//
		resetParentalIndex(newPop, "father_idx", "mother_idx", "oldindex");
		return newPop;
	}

}
