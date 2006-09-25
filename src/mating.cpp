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

#include "mating.h"

namespace simuPOP
{
	mating::mating(double numOffspring,
		PyObject* numOffspringFunc,
		UINT maxNumOffspring,
		UINT mode,
		vectorlu newSubPopSize,
		string newSubPopSizeExpr,
		PyObject* newSubPopSizeFunc)
		:m_numOffspring(numOffspring),
		m_numOffspringFunc(NULL),
		m_maxNumOffspring(maxNumOffspring),
		m_mode(mode), m_firstOffspring(true),
		m_subPopSize(newSubPopSize),
		m_subPopSizeExpr(newSubPopSizeExpr,""),
		m_subPopSizeFunc(NULL)
	{
		DBG_FAILIF( !m_subPopSizeExpr.empty() && newSubPopSizeFunc != NULL,
			ValueError, "Please only specify one of newSubPopSizeExpr and newSubPopSizeFunc.");

		if ( numOffspringFunc != NULL)
		{
			if(!PyCallable_Check(numOffspringFunc))
				throw ValueError("Passed variable is not a callable python function.");

			Py_INCREF(numOffspringFunc);
			m_numOffspringFunc = numOffspringFunc;
		}

		if( newSubPopSizeFunc != NULL)
		{
			if(!PyCallable_Check(newSubPopSizeFunc))
				throw ValueError("Passed variable is not a callable python function.");

			Py_INCREF(newSubPopSizeFunc);
			m_subPopSizeFunc = newSubPopSizeFunc;
		}

		DBG_FAILIF(mode == MATE_BinomialDistribution && maxNumOffspring < 2,
			ValueError, "If mode is MATE_BinomialDistribution, maxNumOffspring should be > 1");
	}


	UINT mating::numOffspring(int gen )
	{
		static double numOS = 0.;

		// use the same numOffspings each generation
		if( ! m_firstOffspring && m_mode == MATE_NumOffspring )
			return static_cast<UINT>(numOS);

		if(m_numOffspringFunc == NULL)
			numOS = m_numOffspring;
		else if( m_mode == MATE_NumOffspring ||	  // case 1, first time
			// case 2, all the time; case 3,4,5, first time
			m_mode == MATE_NumOffspringEachFamily ||
			( m_firstOffspring &&
			( m_mode == MATE_GeometricDistribution  ||
			m_mode == MATE_PoissonDistribution   ||
			m_mode == MATE_BinomialDistribution )))
		{
			PyCallFunc(m_numOffspringFunc, "(i)", gen, numOS, PyObj_As_Double);
		}

		m_firstOffspring = false;

		if( m_mode == MATE_NumOffspring || m_mode == MATE_NumOffspringEachFamily)
		{
			DBG_FAILIF( numOS < 1, ValueError, "Need at least one offspring.");
			return static_cast<UINT>(numOS);
		}
		else if( m_mode == MATE_GeometricDistribution)
		{
			DBG_FAILIF( fcmp_lt(numOS, 0) || fcmp_gt(numOS, 1.), ValueError,
				"P for a geometric distribution should be within [0,1], given " + toStr(numOS));
			UINT nos = rng().randGeometric(numOS);
			return nos;
		}
		else if( m_mode == MATE_PoissonDistribution)
		{
			// FIXME: numOS is no longer the average.
			UINT nos = rng().randPoisson(numOS)+1;
			return nos;
		}
		else if( m_mode == MATE_BinomialDistribution)
		{
			DBG_FAILIF( fcmp_lt(numOS, 0) || fcmp_gt(numOS, 1.), ValueError,
				"P for a Bionomial distribution should be within [0,1], given " + toStr(numOS));
			DBG_FAILIF( m_maxNumOffspring < 1, ValueError,
				"Max number of offspring should be greater than 1. Given "
				+ toStr(m_maxNumOffspring));
			UINT nos = rng().randBinomial(m_maxNumOffspring-1, numOS)+1;
			return nos;
		}
		else
			throw ValueError("Wrong mating numoffspring mode. Should be one of \n"
				"MATE_NumOffspring, MATE_NumOffspringEachFamily and MATE_GEometricDistribution");
	}

	bool mating::formOffGenotype(const vector<Operator* >& ops)
	{
		bool res = true;
		for( vector<Operator *>::const_iterator iop = ops.begin(),
			iopEnd = ops.end(); iop != iopEnd;  ++iop)
		{
			if( (*iop)->formOffGenotype() )
			{
				DBG_DO(DBG_SIMULATOR, cout << "Operator " << (*iop)->__repr__() <<
					" will form offspring" << endl);

				res= false;
				break;
			}
		}
		return res;
	}

	void mating::prepareScratchPop(population& pop, population& scratch)
	{
		/// force new subPop size?
		if( m_subPopSize.empty() && m_subPopSizeExpr.empty()
			&& m_subPopSizeFunc == NULL )
		{
			// use population structure of pop
			// genotype of scratch pop will not be kept.
			scratch.setSubPopStru(pop.subPopSizes(), true);
		}
		else if(! m_subPopSize.empty())			  // set subPoplation size
		{
			// allow  change of pop size of scratch
			scratch.setSubPopStru(this->m_subPopSize, true);

			DBG_FAILIF(scratch.numSubPop() != pop.numSubPop() ,
				ValueError, "number of subPopulaitons must agree. \n Pre: "
				+ toStr( pop.numSubPop()) + " now: " + toStr(scratch.numSubPop() ));
		}
		else if(! m_subPopSizeExpr.empty())		  // evaluate subpopulation size
		{
			m_subPopSizeExpr.setLocalDict( pop.dict());
			vectorf sizef = m_subPopSizeExpr.valueAsArray();
			vectorlu sz(sizef.size());

			for(size_t i=0, iEnd = sizef.size(); i<iEnd; i++)
				sz[i] = static_cast<ULONG>(sizef[i]);

			DBG_DO(DBG_SIMULATOR, cout << "New subpop size " << sz << endl);

			// allow population size change
			scratch.setSubPopStru(sz, true);
			DBG_FAILIF( scratch.numSubPop() != pop.numSubPop() ,
				ValueError, "number of subPopulaitons must agree.\n Pre: "
				+ toStr( pop.numSubPop()) + " now: " + toStr(scratch.numSubPop() ));
		}
		else									  // use m_subPopSizeFunc
		{
			// get generation number
			int gen = pop.gen();
			//          int gen = mainVars().getVarAsInt("gen");
			// convert current pop size to a tuple
			PyObject* curSize = PyTuple_New(pop.numSubPop());

			DBG_ASSERT(curSize!=NULL, SystemError, "Can not convert current pop size to a list");

			for(size_t i=0; i<pop.numSubPop(); ++i)
				PyTuple_SetItem(curSize, i, PyInt_FromLong(pop.subPopSize(i)));

			vectorf res;
			PyCallFunc2(m_subPopSizeFunc, "(iO)", gen, curSize,
				res, PyObj_As_Array);
			Py_XDECREF(curSize);				  // just to be safe

			vectorlu sz(res.size());

			for(size_t i=0; i< res.size(); i++)
				sz[i] = static_cast<ULONG>(res[i]);

			DBG_DO(DBG_SIMULATOR, cout << "New subpop size " << sz << endl);

			// allow  change of pop size of scratch
			scratch.setSubPopStru(sz, true);
			DBG_FAILIF( scratch.numSubPop() != pop.numSubPop() ,
				ValueError, "number of subPopulaitons must agree.\n Pre: "
				+ toStr( pop.numSubPop()) + " now: " + toStr(scratch.numSubPop() ));
		}
	}

	bool noMating::mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit)
	{
		// apply during mating operators
		if( ! ops.empty() )
		{
			for(population::IndIterator it = pop.indBegin(), itEnd = pop.indEnd(); it != itEnd;  ++it)
			{
				for( vector<Operator *>::iterator iop = ops.begin(), iopEnd = ops.end(); iop != iopEnd;  ++iop)
				{
					(*iop)->applyDuringMating(pop, it, NULL, NULL);
				}								  // all during-mating operators
			}
		}
		return true;
	}

	bool binomialSelection::mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit)
	{
		this->resetNumOffspring();
		// scrtach will have the right structure.
		this->prepareScratchPop(pop, scratch);

		DBG_DO(DBG_MATING, m_famSize.clear());

		DBG_ASSERT( pop.numSubPop() == scratch.numSubPop(), SystemError,
			"Number of subpopulation can not be changed.");

		GappedInfoIterator fitness;
        bool selectionOn = pop.hasVar("selection") and pop.getVarAsBool("selection");
        if (selectionOn)
        {
            UINT fit_id = pop.infoIdx("fitness");
	    	fitness = pop.infoBegin(fit_id);
        }

		/// determine if mate() will generate offspring genotype
		bool formOffGeno = this->formOffGenotype(ops);

		// for each subpopulation
		for(UINT sp=0; sp < pop.numSubPop(); ++sp)
		{
			UINT spSize = pop.subPopSize(sp);
			if( spSize == 0 ) continue;

			// if selection is on
			if( selectionOn )
				m_sampler.set( vectorf(fitness+pop.subPopBegin(sp),
					fitness+pop.subPopEnd(sp) ) );

			// choose a parent and genreate m_numOffspring offspring
			ULONG spInd = 0;
			ULONG spIndEnd = scratch.subPopSize(sp);
			while( spInd < spIndEnd)
			{
				individual * parent;
				// choose a parent
				if(selectionOn )
					parent = &pop.ind( m_sampler.get(), sp);
				else
					parent = &pop.ind( rng().randInt(spSize), sp);

				// generate m_numOffspring offspring
				UINT numOS, numOSEnd;
				for(numOS=0, numOSEnd = this->numOffspring(pop.gen() ); numOS < numOSEnd;  numOS++)
				{
					population::IndIterator it = scratch.indBegin(sp) + spInd++;

					if(formOffGeno)
						/// use deep copy!!!!!!!
						it->copyFrom(*parent);

					// apply during mating operators
					for( vector<Operator *>::iterator iop = ops.begin(),
						iopEnd = ops.end(); iop != iopEnd;  ++iop)
					{
						try
						{
							// During mating operator might reject this offspring.
							if(!(*iop)->applyDuringMating(pop, it, parent, NULL))
							{
								spInd --;
								numOS --;
								break;
							}
						}
						catch(...)
						{
							cout << "DuringMating operator " << (*iop)->__repr__() << " throws an exception." << endl << endl;
							throw;
						}
					}							  // all during-mating operators
					// success
					if( spInd == spIndEnd )
					{
						numOS++;
						break;
					}
				}								  // offsrping for each parent
				// record family size
				DBG_DO(DBG_MATING, m_famSize.push_back( numOS ));
			}									  // all offspring
		}										  // all subpopulation.

		if(submit)
			submitScratch(pop, scratch);
		return true;
	}

	bool randomMating::mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit)
	{
		bool hasSexChrom = pop.sexChrom();

		this->resetNumOffspring();
		// scrtach will have the right structure.
		this->prepareScratchPop(pop, scratch);

		DBG_DO(DBG_MATING, m_famSize.clear());

		DBG_ASSERT( pop.numSubPop() == scratch.numSubPop(), SystemError,
			"Number of subpopulation can not be changed.");

		// empty fitness means no selection
		GappedInfoIterator fitness;
        bool selectionOn = pop.hasVar("selection") and pop.getVarAsBool("selection");
        if (selectionOn)
        {
            UINT fit_id = pop.infoIdx("fitness");
	    	fitness = pop.infoBegin(fit_id);
        }

		/// determine if any during-mating operator will generate offspring genotype
		/// so mating scheme does not have to do it.
		bool formOffGeno = this->formOffGenotype(ops);

		UINT numMale, numFemale;

		// cache chromBegin, chromEnd for better performance.
		const vectoru & chIdx = pop.chromIndex();

		/// random mating happens within each subpopulation
		for(UINT sp=0; sp < pop.numSubPop(); ++sp)
		{
			ULONG spSize = pop.subPopSize(sp);
			if( spSize == 0 )
				continue;

			numMale = 0;
			for(population::IndIterator it=pop.indBegin(sp), itEnd = pop.indEnd(sp); it < itEnd; ++it)
				if(it->sex() == Male)
					numMale ++;

			// allocate memory at first for performance reasons
			m_maleIndex.resize(numMale);
			m_femaleIndex.resize(spSize-numMale);

			numMale = 0;
			numFemale = 0;

			for( ULONG it = pop.subPopBegin(sp), itEnd = pop.subPopEnd(sp); it < itEnd;  it++)
			{
				if( pop.ind(it).sex() == Male)
					m_maleIndex[numMale++] = it;
				else
					m_femaleIndex[numFemale++] = it;
			}

			/// now, all individuals of needToFind sex is collected
			if( (numMale == 0 || numFemale ==0 ) && !m_contWhenUniSex )
			{
				throw ValueError("Subpopulation becomes uni-sex. Can not continue. \n"
					"You can use ignoreParentsSex (do not check parents' sex) or \ncontWhenUnixSex "
					"(same sex mating if have to) options to get around this problem.");
			}

			DBG_ASSERT( numFemale + numMale == spSize, SystemError,
				"Wrong number of male/female.");

			/// if selection is on
			if( selectionOn)
			{
				m_maleFitness.resize(numMale);
				m_femaleFitness.resize(numFemale);

				size_t ind;

				for( ind = 0; ind < numMale; ++ind)
					m_maleFitness[ind] = fitness[ m_maleIndex[ind] ];
				for( ind = 0; ind < numFemale; ++ind)
					m_femaleFitness[ind] = fitness[ m_femaleIndex[ind] ];

				m_malesampler.set(m_maleFitness);
				m_femalesampler.set(m_femaleFitness);
			}

			// generate scratch.subPopSize(sp) individuals.
			ULONG spInd = 0;
			ULONG spIndEnd = scratch.subPopSize(sp);

			// in case that mating will form genotype and we need
			// free recombination between chromosomes.
			BernulliTrials bt(rng(), vectorf(1, 0.5), 2*pop.numChrom());

			while( spInd < spIndEnd)
			{
				// randomly choose parents
				individual * dad, *mom;
				RNG& rnd = rng();

				if( selectionOn)			  // with selection
				{
					// using weidhted sampler.
					if( numMale != 0 )
						dad = &pop.ind( m_maleIndex[ m_malesampler.get() ] );
					else
						dad = &pop.ind( m_femaleIndex[ m_femalesampler.get() ] );

					if( numFemale != 0 )
						mom = &pop.ind( m_femaleIndex[ m_femalesampler.get() ] );
					else
						mom = &pop.ind( m_maleIndex[ m_malesampler.get() ]);
				}
				else
				{
					// using random sample.
					if( numMale != 0 )
						dad = &pop.ind( m_maleIndex[ rnd.randInt(numMale) ]);
					else
						dad = &pop.ind( m_femaleIndex[ rnd.randInt(numFemale) ]);

					if( numFemale != 0 )
						mom = &pop.ind( m_femaleIndex[ rnd.randInt(numFemale) ]);
					else
						mom = &pop.ind( m_maleIndex[ rnd.randInt(numMale) ]);
				}

				// generate m_numOffspring offspring per mating
				UINT numOS=0, numOSEnd;
				for(numOS=0, numOSEnd = this->numOffspring(pop.gen()); numOS < numOSEnd;  numOS++)
				{
					population::IndIterator it = scratch.indBegin(sp) + spInd++;

					//
					// assign sex randomly
					if( ! hasSexChrom)
					{
						int offSex = rnd.randInt(2);
						it->setSex( offSex==0?Male:Female);
					}

					if( formOffGeno )			  // use the default no recombination random mating.
					{
						bt.trial();
						const BitSet& bs = bt.succ(0);

						// initialize to avoid compiler complains
						int dadPloidy=0, momPloidy=0;
						GenoIterator cd[2], cm[2], offd, offm;
						cd[0] = dad->genoBegin(0);
						cd[1] = dad->genoBegin(1);
						cm[0] = mom->genoBegin(0);
						cm[1] = mom->genoBegin(1);
						offd = it->genoBegin(0);
						offm = it->genoBegin(1);

						for(UINT ch=0, chEnd = dad->numChrom(); ch < chEnd; ++ch)
						{
							// bs is 2*totNumLoci() long
							dadPloidy = bs[ch];
							momPloidy = bs[ch+chEnd];
							for(size_t gt = chIdx[ch]; gt < chIdx[ch+1]; ++gt)
							{
								offd[gt] = cd[dadPloidy][gt];
								offm[gt] = cm[momPloidy][gt];
							}
						}

						// last chromosome (sex chromosomes) determine sex
						if( hasSexChrom )
						{
							if( dadPloidy == 1)	  // Y of XY
								it->setSex(Male);
							else
								it->setSex(Female);
						}
					}

					/// apply all during mating operators
					for( vector<Operator *>::iterator iop = ops.begin(), iopEnd = ops.end(); iop != iopEnd;  ++iop)
					{
						try
						{
							// During mating operator might reject this offspring.
							if(!(*iop)->applyDuringMating(pop, it, dad, mom))
							{
								spInd --;
								numOS --;
								break;
							}
						}
						catch(...)
						{
							cout << "DuringMating operator " << (*iop)->__repr__() << " throws an exception." << endl << endl;
							throw;
						}
					}
					// appy operators
					// success
					if( spInd == spIndEnd )
					{
						numOS++;
						break;
					}
				}								  // each offspring
				// and then break here since spInd == spIndEnd.
				// record family size
				DBG_DO(DBG_MATING, m_famSize.push_back( numOS ));
			}
		}										  // each subPop

		if(submit)
			submitScratch(pop, scratch);
		return true;
	}

	vectorf FreqTrajectoryStoch(ULONG curGen, double freq, long N,
		PyObject* NtFunc, vectorf fitness, PyObject* fitnessFunc,
		ULONG minMutAge, ULONG maxMutAge, int ploidy,
        bool restartIfFail, long maxAttempts, bool allowFixation)
	{
		if( curGen >0 && minMutAge > curGen )
			minMutAge = curGen;
		if( curGen >0 && maxMutAge > curGen )
			maxMutAge = curGen;

		DBG_FAILIF( maxMutAge < minMutAge, ValueError, "maxMutAge should >= minMutAge");
		DBG_FAILIF( curGen == 0 && (NtFunc != NULL || fitnessFunc != NULL),
			ValueError, "curGen should be > 0 if NtFunc or fitnessFunc is defined.");
		DBG_FAILIF( curGen > 0 && curGen < maxMutAge, ValueError,
			"curGen should be >= maxMutAge");

		// is NtFunc callable?
		if( NtFunc != NULL )
		{
			if( ! PyCallable_Check(NtFunc) )
				throw ValueError("NtFunc is not a valid Python function.");
			else
				// increase the ref, just to be safe
				Py_INCREF(NtFunc);
		}

		// 1, 1+s1, 1+s2
		double s1, s2;
		if( fitnessFunc!= NULL )
		{
			if( ! PyCallable_Check(fitnessFunc) )
				throw ValueError("fitnessFunc is not a valid Python function.");
			else
				// increase the ref, just to be safe
				Py_INCREF(fitnessFunc);
		}
		else if( fitness.empty() )
			// default to [1,1,1]
		{
			s1 = 0.;
			s2 = 0.;
		}
		else if( fitness.size() != 3 || fitness[0] == 0)
		{
			throw ValueError("s should be a vector of length 3. (for AA, Aa and aa)");
		}
		else
		{
			// convert to the form 1, s1, s2
			s1 = fitness[1] / fitness[0] - 1.;
			s2 = fitness[2] / fitness[0] - 1.;
		}

		// get current population size
		vectori Ntmp(1, N);
		if( NtFunc != NULL)
		{
			PyCallFunc(NtFunc, "(i)", curGen, Ntmp, PyObj_As_IntArray);
			DBG_ASSERT( Ntmp.size() >= 1, ValueError,
				"Return value from NtFunc should be an array of size >= 1");
			// Ntmp[0] will be the total size.
			for(size_t i=1; i<Ntmp.size(); ++i)
				Ntmp[0] += Ntmp[i];
		}

		// all calculated population size
		vectorlu Nt(1, Ntmp[0]);
		// copies of allele a at each genertion.
		vectorlu it(1, static_cast<long>(ploidy*Ntmp[0]*freq));
		// allele frequency of allele a at each geneation
		vectorf xt(1, freq);
		// store calculated fitness s1, s2, if necessary
		// note that s_T will never be used.
		vectorf s1_cache(1,0), s2_cache(1,0);
		// will be used to store return value of sFunc
		vectorf s_vec;

		// t is current generation number.
		ULONG idx = 0;

		// a,b,c etc for solving the quadratic equation
		long failedCount = 0;
		long tooLongCount = 0;
		long tooShortCount = 0;
		long invalidCount = 0;
		// a,b,c etc for solving the quadratic equation
		double a,b,c,b2_4ac,y1,y2,y;
		while( true )
		{
			// if too many fails
			if( failedCount >= maxAttempts)
				break;

			// first get N(t-1), if it has not been calculated
			if( idx+1 >= Nt.size() )
			{
				if( NtFunc != NULL)
				{
					PyCallFunc(NtFunc, "(i)", curGen-idx-1, Ntmp, PyObj_As_IntArray);
					DBG_ASSERT( Ntmp.size() >= 1, ValueError,
						"Return value from NtFunc should be an array of size >= 1");
					// Ntmp[0] will be the total size.
					for(size_t i=1; i<Ntmp.size(); ++i)
						Ntmp[0] += Ntmp[i];
				}
				Nt.push_back(Ntmp[0]);
			}
			//
			// get fitness
			if( fitnessFunc != NULL )
			{
				if( idx+1 >= s1_cache.size() )
				{
					PyCallFunc(fitnessFunc, "(i)", curGen-idx-1, s_vec, PyObj_As_Array);

					DBG_ASSERT(s_vec.size()==3 || s_vec[0] != 0., ValueError,
						"Returned value from sFunc should be a vector of size 3");

					s1 = s_vec[1]/s_vec[0] - 1.;
					s2 = s_vec[2]/s_vec[0] - 1.;
					s1_cache.push_back(s1);
					s2_cache.push_back(s2);
				}
				else
				{
					s1 = s1_cache[idx];
					s2 = s2_cache[idx];
				}
			}
			// given x(t)
			// calculate y=x(t-1)' by solving an equation
			//
			//  x_t = y(1+s2 y+s1 (1-y))/(1+s2 y+2 s1 y(1-y))
			//
			if( s1 == 0. && s2 == 0. )
			{
				// special case when a = 0.
				y = xt[idx];
			}
			else
			{
				a = s2*xt[idx] - 2*s1*xt[idx] - s2 + s1;
				b = 2*s1*xt[idx] - 1 - s1;
				c = xt[idx];
				b2_4ac = b*b-4*a*c;

				if( b2_4ac < 0 )
					throw ValueError("Quadratic function does not yield a valid solution");

				y1 = (-b+sqrt(b2_4ac))/(2*a);
				y2 = (-b-sqrt(b2_4ac))/(2*a);

				// choose one of the solutions
				if( y1 < 0. || y1 > 1.0 )
				{
					if( y2 >= 0. && y2 <= 1.0)
						y = y2;
					else
						throw ValueError("None of the solutions is valid");
				}
				else
				{
					// if y1 is valid
					if( y2 >= 0. && y2 <= 1.0)
						throw ValueError("Both solutions are valid. Further decision is needed. y1: " +
							toStr(y1) + " y2: " + toStr(y2) + " xt: " + toStr(xt[idx]));
					else
						y = y1;
				}
			}

			if( it.size() < idx+2 )
			{
				it.push_back(0);
				xt.push_back(0);
			}
			// y is obtained, is the expected allele frequency for the next generation t-1
			it[idx+1] = rng().randBinomial(ploidy*Nt[idx+1], y);
			xt[idx+1] = it[idx+1]/static_cast<double>(ploidy*Nt[idx+1]);

			if( it[idx+1] == 0 )
			{
				if( idx < minMutAge)
				{
					// cout << "Path too short. Retrying" << endl;
					failedCount ++;
					tooShortCount ++;
					idx = 0;
					continue;
				}
				// need 0, 1, ....,
				else if( it[idx] == 1 || allowFixation)
				{
					// success
					break;
				}
				else
				{
					DBG_DO(DBG_MATING, cout << "Reaching 0, but next gen has more than 1 allele a" << endl);
					invalidCount ++;
					failedCount ++;
					// restart
					idx = 0;
				}
			}
			else if( it[idx+1] == ploidy*Nt[idx+1] )
				// when the allele get fixed, restart
			{
				if(allowFixation)
				{
					break;
				} else {
					failedCount ++;
					invalidCount ++;
					idx = 0;
				}
			}
			// if not done, but t already reaches T
			else if( idx == maxMutAge )
			{
				if( restartIfFail )
				{
					idx = 0;
					failedCount ++;
					tooLongCount ++;
					DBG_DO(DBG_GENERAL, cout << "Warning: reaching max gnerations " + toStr(maxMutAge) + ". Restart the process." << endl);
					continue;
				}
				else
				{
					//DBG_DO(DBG_GENERAL, cout << "Warning: reaching max gnerations. Return whatever I have now." << endl);
					break;
				}
			}
			else
				// go to next generation
				idx++;
		}
		// report potential problems
		if( tooLongCount > 0 )
			cout << "Trajectories regenerated due to long (> " << maxMutAge << ") path: " << tooLongCount << " times." << endl;
		if( tooShortCount > 0 )
			cout << "Trajectories regenerated due to short (< " << minMutAge << ") path: " << tooShortCount << " times." << endl;
		if( invalidCount > 0 )
			cout << "Trajectories regenerated due to invalid path: " << invalidCount << " times. " << endl;

		// clean up
		if( NtFunc != NULL)
			Py_DECREF(NtFunc);
		if( fitnessFunc != NULL)
			Py_DECREF(fitnessFunc);

		// number of valid generation is idx+1
		vectorf traj(idx+1);
		if( failedCount < maxAttempts )
		{
			for(ULONG i=0; i<=idx; ++i)
				traj[i] = xt[idx-i];
		}
		return traj;
	}


	matrix FreqTrajectoryMultiStoch( ULONG curGen,
		vectorf freq, long N,
		PyObject* NtFunc, vectorf fitness, PyObject* fitnessFunc,
		ULONG minMutAge, ULONG maxMutAge, int ploidy,
        bool restartIfFail, long maxAttempts)
	{
		size_t nLoci = freq.size();
		size_t i, j, curI, nextI;

		if( curGen >0 && minMutAge > curGen )
			minMutAge = curGen;
		if( curGen >0 && maxMutAge > curGen )
			maxMutAge = curGen;

		DBG_ASSERT( minMutAge <= maxMutAge, ValueError, "minMutAge should be <= maxMutAge. ");
		DBG_ASSERT( nLoci > 0, ValueError, "Number of loci should be at least one");
		DBG_FAILIF( curGen == 0 && (NtFunc != NULL || fitnessFunc != NULL),
			ValueError, "curGen should be > 0 if NtFunc or fitnessFunc is defined.");
		DBG_FAILIF( curGen > 0 && curGen < maxMutAge, ValueError,
			"curGen should be >= maxMutAge");

		matrix result(nLoci);

		// in the cases of independent and constant selection pressure
		// easy case.
		if( fitnessFunc == NULL)
		{
			DBG_FAILIF( (!fitness.empty()) && (fitness.size() != nLoci*3),
				ValueError, "Wrong s length " + toStr(fitness.size()));

			for( i=0; i<nLoci; ++i)
			{
				if( ! fitness.empty() )
					result[i] = FreqTrajectoryStoch(curGen, freq[i], N, NtFunc,
						vectorf(fitness.begin()+3*i, fitness.begin()+3*(i+1)),
						NULL, minMutAge, maxMutAge, ploidy, restartIfFail, maxAttempts);
				else
					result[i] = FreqTrajectoryStoch(curGen, freq[i], N, NtFunc,
						vectorf(), NULL, minMutAge, maxMutAge, ploidy, restartIfFail, maxAttempts);
			}
			return result;
		}

		// other wise, use a vectorized version of ...
		// is NtFunc callable?
		if( NtFunc != NULL )
		{
			if( !PyCallable_Check(NtFunc) )
				throw ValueError("NtFunc is not a valid Python function.");
			else
				// increase the ref, just to be safe
				Py_INCREF(NtFunc);
		}

		// sAll will store returned value of sFunc,
		// and be converted to 1, 1+s1, 1+s2 format ...
		vectorf sAll;
		if( fitnessFunc != NULL )
		{
			if( !PyCallable_Check(fitnessFunc) )
				throw ValueError("sFunc is not a valid Python function.");
			else
				// increase the ref, just to be safe
				Py_INCREF(fitnessFunc);
		}
		else if( fitness.empty() )
			// default to [1,1,1]
		{
			sAll.resize(nLoci*3, 0.);
		}
		else if( fitness.size() != 3*nLoci)
		{
			throw ValueError("s should be a vector of length 3 times number of loci. (for AA, Aa and aa etc)");
		}
		else
		{
			for(i=0; i<nLoci; ++i)
			{
				// convert to the form 1, s1, s2
				sAll[3*i+1] = fitness[3*i+1] / fitness[3*i] - 1.;
				sAll[3*i+2] = fitness[3*i+2] / fitness[3*i] - 1.;
				sAll[3*i] = 0.;
			}
		}

		// get current population size
		vectori Ntmp(1, N);
		if( NtFunc != NULL)
		{
			PyCallFunc(NtFunc, "(i)", curGen, Ntmp, PyObj_As_IntArray);
			DBG_ASSERT( Ntmp.size() >= 1, ValueError,
				"Return value from NtFunc should be an array of size >= 1");
			// Ntmp[0] will be the total size.
			for(size_t i=1; i<Ntmp.size(); ++i)
				Ntmp[0] += Ntmp[i];
		}

		// all calculated population size
		vectorlu Nt(1, Ntmp[0]);
		// copies of allele a at each genertion.
		vectorlu it(nLoci);
		for(i=0; i<nLoci; ++i)
			it[i] = static_cast<long>(ploidy*Ntmp[0]*freq[i]);
		// allele frequency of allele a at each geneation
		vectorf xt(nLoci);
		for(i=0; i<nLoci; ++i)
			xt[i] = freq[i];

		ULONG idx = 0;

		long failedCount = 0;
		long tooLongCount = 0;
		long tooShortCount = 0;
		long invalidCount = 0;
		// a,b,c etc for solving the quadratic equation
		double s1,s2,x,a,b,c,b2_4ac,y1,y2,y;
		// whether or not each locus is done
		vector<bool> done(nLoci, false);
		//
		while( true )
		{
			// if too many fails
			if( failedCount >= maxAttempts)
				break;

			// first get N(t-1), if it has not been calculated
			if( idx+1 >= Nt.size() )
			{
				if( NtFunc != NULL)
				{
					PyCallFunc(NtFunc, "(i)", curGen-idx-1, Ntmp, PyObj_As_IntArray);
					DBG_ASSERT( Ntmp.size() >= 1, ValueError,
						"Return value from NtFunc should be an array of size >= 1");
					// Ntmp[0] will be the total size.
					for(size_t i=1; i<Ntmp.size(); ++i)
						Ntmp[0] += Ntmp[i];
				}
				Nt.push_back(Ntmp[0]);
			}
			//
			// get fitness, since it will change according to
			// xt, I do not cache the result
			if( fitnessFunc != NULL )
			{
				// compile allele frequency... and pass
				PyObject* freqObj = Double_Vec_As_NumArray( xt.begin()+nLoci*idx, xt.begin()+nLoci*(idx+1) );
				PyCallFunc2(fitnessFunc, "(iO)", curGen-idx-1, freqObj, sAll, PyObj_As_Array);

				DBG_ASSERT(sAll.size()==3*nLoci, ValueError,
					"Returned value from sFunc should be a vector of size 3");

				for(i=0; i<nLoci; ++i)
				{
					// convert to the form 1, s1, s2
					sAll[3*i+1] = sAll[3*i+1] / sAll[3*i] - 1.;
					sAll[3*i+2] = sAll[3*i+2] / sAll[3*i] - 1.;
					sAll[3*i] = 0.;
				}
			}

			bool restart = false;
			// handle each locus
			for(i = 0; i<nLoci; ++i)
			{
				curI = idx*nLoci+i;
				nextI = (idx+1)*nLoci+i;

				// allocate space
				if( it.size() < (idx+2)*nLoci )
				{
					for(j=0;j<nLoci;++j)
					{
						it.push_back(0);
						xt.push_back(0.);
					}
				}
				// if done
				if( done[i] )
				{
					it[nextI] = 0;
					xt[nextI] = 0.;
					continue;
				}
				// given x(t)
				// calculate y=x(t-1)' by solving an equation
				//
				//  x_t = y(1+s2 y+s1 (1-y))/(1+s2 y+2 s1 y(1-y))
				//
				s1 = sAll[i*3+1];
				s2 = sAll[i*3+2];
				x  = xt[curI];
				if( s1 == 0. && s2 == 0. )
				{
					// special case when a = 0.
					y = x;
				}
				else
				{
					a = s2*x - 2*s1*x - s2 + s1;
					b = 2*s1*x - 1 - s1;
					c = x;
					b2_4ac = b*b-4*a*c;

					if( b2_4ac < 0 )
						throw ValueError("Quadratic function does not yield a valid solution");

					y1 = (-b+sqrt(b2_4ac))/(2*a);
					y2 = (-b-sqrt(b2_4ac))/(2*a);

					// choose one of the solutions
					if( y1 < 0. || y1 > 1.0 )
					{
						if( y2 >= 0. && y2 <= 1.0)
							y = y2;
						else
							throw ValueError("None of the solutions is valid");
					}
					else
					{
						// if y1 is valid
						if( y2 >= 0. && y2 <= 1.0)
							throw ValueError("Both solutions are valid. Further decision is needed. y1: " +
								toStr(y1) + " y2: " + toStr(y2) + " xt: " + toStr(xt[idx]));
						else
							y = y1;
					}
				}

				// y is obtained, is the expected allele frequency for the next generation t-1
				it[nextI] = rng().randBinomial(ploidy*Nt[idx+1], y);
				xt[nextI] = it[nextI]/static_cast<double>(ploidy*Nt[idx+1]);

				if( it[nextI] == 0 )
				{
					if( idx < minMutAge )
					{
						tooShortCount ++;
						restart = true;
						break;
					}
					// need 0, 1, ...., good...
					if( it[curI] == 1 )
					{
						done[i] = true;
					}
					else
					{
						// restart
						invalidCount ++;
						restart = true;
						break;
					}
				}
				else if( it[nextI] == ploidy*Nt[idx+1] )
				{
					// when the allele get fixed, restart
					invalidCount ++;
					restart = true;
					break;
				}
			}									  // end of for each locus

			// break from inside
			if( restart || (idx==maxMutAge && restartIfFail))
			{
				failedCount ++;
				idx = 0;
				for( j=0; j<nLoci; ++j)
					done[j] = false;
				if(idx == maxMutAge)
					tooLongCount ++;
				continue;
			}

			// if all done?
			if( find(done.begin(), done.end(), false) == done.end() )
				break;

			DBG_DO(DBG_DEVEL, cout << idx << " freq= " <<
				vectorf(xt.begin()+idx*nLoci, xt.begin()+(idx+1)*nLoci) <<
				" s= " << sAll << endl);
			//
			// if not done, but t already reaches T
			if( idx == maxMutAge )
			{
				DBG_DO(DBG_MATING, cout << "Warning: reaching T generations. Return whatever I have now." << endl);
				break;
			}

			// go to next generation
			idx ++;
		}
		// report potential problems
		if( tooLongCount > 0 )
			cout << "Trajectories regenerated due to long path: " << tooLongCount << " times. " << endl;
		if( tooShortCount > 0 )
			cout << "Trajectories regenerated due to short path: " << tooShortCount << " times. " << endl;
		if( invalidCount > 0 )
			cout << "Trajectories regenerated due to invalid path: " << invalidCount << " times. " << endl;

		// clean up
		if( NtFunc != NULL)
			Py_DECREF(NtFunc);
		if( fitnessFunc != NULL)
			Py_DECREF(fitnessFunc);

		// number of valid generation is idx+1
		if( failedCount < maxAttempts )
		{
			vectorf traj(idx+1);
			for(i=0; i<nLoci; ++i)
			{
				for(j=0; j<=idx; ++j)
					traj[j] = xt[nLoci*(idx-j)+i];
				for(j=0; j<traj.size() && traj[j]==0.; j++);
				result[i] = vectorf(traj.begin()+j, traj.end());
			}
		}
		return result;
	}


	// simulate trajectory
	vectorf FreqTrajectorySelSim(
		double sel,								  // strength of selection coef  ::8
		long Ne,								  // effective population size ::9
		double freq,							  // initial freq ::10
		double dom_h,							  // strength of dominance ::27
		int selection							  // selection ::5
		)
	{
		//
		// first step, prepeare u:
		//
		// probability of zero for traj condition
		// originally prep_table
		long   N = Ne;
		double s = sel*4.0;
		double h = dom_h;

		vector<double> u(N+1);
		double ratio = ((1-s/4.0)/(1+s/4.0));
		vector<double> prod(N);
		prod[0] = 1;

		for(long i=1;i<N;i++)
		{
			if(selection==1)
				prod[i] = prod[i-1]*ratio;
			if(selection==2)
				prod[i] = prod[i-1]*((1-(s)*(1-2*((double)i/N)))/(1+(s)*(1-2*((double)i/N))));
			if(selection==3)
				prod[i] = prod[i-1]*( 1 - (s/2.0)*(((double)i/N)+h*(1-2*((double)i/N)) ) )/
					( 1 + (s/2.0)*(((double)i/N)+h*(1-2*((double)i/N)) ) );
		}
		u[N-1] = prod[N-1];
		for(long i = N-2; i>0; i--)
			u[i] = u[i+1] + prod[i];
		double first= u[1];
		for(long i=1; i<N; i++)
			u[i] /= (1.0 + first);

		//
		// second step
		//
		// prepare log look up table, necessary?
		vector<double> log_lookup(N+1);
		for(int i=1;i<=N;i++)
			log_lookup[i] = -(double)N*log((double)i/(i+1));

		//
		// third step
		//
		// start
		int jump=-1;
		double weight_up;

		int j = int(double(N)*freq);
		if((j == 0)||(j>=N))
			throw " Starting number selected types = 0 or >=N ";

		double * pos_log_ptr = &log_lookup[j];
		double * neg_log_ptr = &log_lookup[N-j-1];
		double * u0_pnt_j   = &u[j];
		double * u0_pnt_jp1 = &u[j+1];

		double lambda_driving = 1+s/4.0;
		double tot_birth = 2.0;
		double prob_up = lambda_driving/tot_birth;

		// accumulative t_time
		double past_traj_time = 0;

		vector<double> traj_time(1, 0);
		vector<double> traj_freq(1, double(j)/N);
		vector<double> traj_integral_pos(1, 0);
		vector<double> traj_integral_neg(1, 0);

		double prev_pos  = 0;
		double prev_neg  = 0;
		double t_time = 0;

		traj_time.push_back(0);
		traj_freq.push_back(0);

		double position_coeff;

		while(true)
		{
			traj_integral_pos.back() = prev_pos + t_time*(*pos_log_ptr);
			traj_integral_neg.back() = prev_neg + t_time*(*neg_log_ptr);

			position_coeff = (double)j * (N-j)/N;

			//Removed Factor of 2
			t_time = -(1/position_coeff)*log(rng().randUniform01() ) /double(N);
			past_traj_time += t_time;
			weight_up = (*u0_pnt_jp1)/(*u0_pnt_j);

			if(selection==2)
				prob_up = (1+(s)*(1-2*((double)j/N))) / 2.0;
			if(selection==3)
				prob_up = ( 1 + (s/2.0)*(((double)j/N) +
					h*(1-2*(double(j)/N)) ) ) / 2.0 ;

			//Jump
			if( rng().randUniform01() < prob_up*weight_up)
			{
				if(jump==1)
				{
					pos_log_ptr++;
					neg_log_ptr--;
				}
				j++;
				u0_pnt_jp1++;
				u0_pnt_j++;
				jump = 1;
			}
			else
			{
				if(jump==-1)
				{
					pos_log_ptr--;
					neg_log_ptr++;
				}
				j--;
				u0_pnt_jp1--;
				u0_pnt_j--;
				jump=-1;
			}

			prev_pos = traj_integral_pos.back();
			prev_neg = traj_integral_neg.back();
			traj_time.back() = past_traj_time;
			traj_freq.back() = (double)j/N;

			if ( j == 0)
				break;
			if ( j == N )
				cout << "Reach allele N, this should be rare" << endl;

			if ( traj_freq.size() % 10000 == 1)
				cout << "Size " << traj_freq.size() << " at freq " << traj_freq.back() << endl;

			traj_freq.push_back(0);
			traj_time.push_back(0);
			traj_integral_pos.push_back(0);
			traj_integral_neg.push_back(0);
		}

		// now we need to translate 4Nt to generations
		// g = 4Nt, g is generation
		double maxTime = traj_time.back();
		vectorf gen_freq( static_cast<size_t>(maxTime * 4 * Ne)+1);
		size_t ngen = gen_freq.size();
		int curTime = traj_time.size()-1;

		for(size_t g=0; g<ngen; ++g)
		{
			double t = (static_cast<size_t>(maxTime*4*Ne) - g)/(4.*Ne);
			// find, from the right, most close to t
			// idea case: x< t <cur
			while( curTime >= 1 && traj_time[curTime-1] > t)
				--curTime;
			gen_freq[g] = traj_freq[curTime-1] +
				(traj_freq[curTime]-traj_freq[curTime-1])*(t-traj_time[curTime-1]);
		}

		return gen_freq;
	}

	vectorf FreqTrajectoryForward(double lowbound, double highbound,
		int disAge, double grate, long N0, double seleCo)
	{
		vector<double> DissamplePath(disAge+1);

		int ftime = disAge;
		// DissamplePath is the current allele frequency
		DissamplePath[0] = 0;

		int trying=0;
		while(DissamplePath[0] <= lowbound || DissamplePath[0] >= highbound)
		{
			if((++trying)%1000==0)
				cout<<"Trying path "<<trying<<" times\n";
			double num = 1;
			//initially 1 copy of mutant allele
			double Nt = N0*exp(-ftime*grate);
			double fre = 1/(2*Nt);

			// initial allele frequency
			DissamplePath[ftime] = fre;

			// backward in array, but forward in time.
			for(int gth = ftime-1; gth>= 0; gth--)
			{
				// population size at generation gth.
				Nt = N0*exp(-grate*gth);
				// number of alleles is less than four, use
				// poisson approximation.
				if(num<4)
				{
					double lamda;
					// expected allele frequency: sp(1-p)
					lamda = (fre + seleCo*fre*(1-fre))*2*Nt;
					num = rng().randPoisson(lamda);
					// restart if num<=0
					if(num<= 0)
					{
						DissamplePath[0] = 0;
						break;
					}
					fre = num/(2*Nt);
					// restart?
					if(fre>= 1)
					{
						DissamplePath[0] = 1;
						break;
					}
					DissamplePath[gth] = fre;
				}
				// diffusion process approximation
				else
				{
					// use uniform approximation???
					double min = seleCo*fre*(1-fre)-sqrt(3*fre*(1-fre)/(2*Nt));
					double max = seleCo*fre*(1-fre)+sqrt(3*fre*(1-fre)/(2*Nt));
					double psv = min + (max-min)*rng().randUniform01();
					fre = fre+psv;
					if(fre<= 0)
					{
						DissamplePath[0] = 0;
						break;
					}
					if(fre>= 1)
					{
						DissamplePath[0] = 1;
						break;
					}
					num = fre* 2* Nt;
					DissamplePath[gth] = fre;
				}
			}
		}										  // while
		// reverse the result and return
		vectorf gen_freq(disAge+1);
		for(int i=0; i<disAge+1; ++i)
			gen_freq[i] = DissamplePath[disAge-i];

		return gen_freq;
	}

	vectorlu controlledMating::countAlleles(population& pop, const vectori& loci, const vectori& alleles)
	{
		vectorlu alleleNum(loci.size(), 0L);
		for(size_t l=0; l < loci.size(); ++l)
		{
			int loc = loci[l];
			Allele ale = alleles[l];
			// go through all alleles
			for(GappedAlleleIterator a=pop.alleleBegin(loc),
				aEnd=pop.alleleEnd(loc); a != aEnd; ++a)
			{
				if( AlleleUnsigned(*a) == ale )
					alleleNum[l]++;
			}
		}
		return alleleNum;
	}

	bool controlledMating::mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit)
	{
		// first call the function and get the range
		vectorf freqRange;

		PyCallFunc( m_freqFunc, "(i)", pop.gen(),
			freqRange, PyObj_As_Array);

		DBG_ASSERT( freqRange.size() == m_loci.size() || freqRange.size() == 2 * m_loci.size(),
			ValueError, "Length of returned range should have the same or double the number of loci.");

		vectorlu alleleNum;

#ifndef OPTIMIZED
		// calculate allele frequen at these loci
		alleleNum = countAlleles(pop, m_loci, m_alleles);

		DBG_DO(DBG_MATING, cout << "Range of allele frequencies at generation "
			<< pop.gen() << " is " << freqRange << endl);
#endif

		// compare integer is easier and more acurate, and we need to make sure
		// there is one allele when the frequency is greater than 0.
		vectorlu alleleRange(2*m_loci.size(), 0);
		if( freqRange.size() == m_loci.size() )
		{
			for(size_t i=0; i<m_loci.size(); ++i)
			{
				if (freqRange[i] < 0. )
					freqRange[i] = 0.;
				if (freqRange[i] > 1. )
					freqRange[i] = 1.;

				alleleRange[2*i] = static_cast<ULONG>(freqRange[i]*pop.popSize()*pop.ploidy());
				if( freqRange[i]>0 && alleleRange[2*i] == 0)
					alleleRange[2*i] = 1;
				// upper bound using range.
				alleleRange[2*i+1] = static_cast<ULONG>((freqRange[i]+m_range)*pop.popSize()*pop.ploidy())+1;
#ifndef OPTIMIZED
				if( alleleNum[i] == 0 && alleleRange[2*i] > 0 )
					throw ValueError("No allele exists so there is no way to reach specified allele frequency.\n"
						"Locus " + toStr(m_loci[i]) + " at generation " + toStr(pop.gen()) );
#endif
			}
		}
		else
			// returned values are [low1, high1, low2, high2...]
		{
			for(size_t i=0; i<m_loci.size(); ++i)
			{
				if (freqRange[2*i] < 0. )
					freqRange[2*i] = 0.;
				if (freqRange[2*i+1] > 1. )
					freqRange[2*i+1] = 1.;

				alleleRange[2*i] = static_cast<ULONG>(freqRange[2*i]*pop.popSize()*pop.ploidy());

				DBG_FAILIF( freqRange[2*i] > freqRange[2*i+1], ValueError,
					"Incorrect frequency range: " + toStr(freqRange[2*i]) +
					" - " + toStr(freqRange[2*i+1]));

				if( freqRange[2*i]>0 && alleleRange[2*i] == 0)
					alleleRange[2*i] = 1;
				alleleRange[2*i+1] = static_cast<ULONG>(freqRange[2*i+1]*pop.popSize()*pop.ploidy())+1;
#ifndef OPTIMIZED
				if( alleleNum[i] == 0 && alleleRange[2*i] > 0 )
					throw ValueError("No allele exists so there is no way to reach specified allele frequency.\n"
						"Locus " + toStr(m_loci[i]) + " at generation " + toStr(pop.gen()) );
#endif
			}
		}

		while(true)
		{
			// do the mating
			m_matingScheme->mate(pop, scratch, ops, false);

			// check allele frequency
			alleleNum = countAlleles(scratch, m_loci, m_alleles);

			DBG_DO(DBG_MATING, cout << "mating finished, new count "
				<< alleleNum << " range " << alleleRange << endl);
			//
			bool succ = true;
			for(size_t i=0; i<m_loci.size(); ++i)
			{
				if( alleleNum[i] < alleleRange[2*i] || alleleNum[i] > alleleRange[2*i+1] )
				{
					succ = false;
					break;
				}
			}
			// cout << "succ " << succ << "submit " << submit << endl;
			if( succ && submit)
			{
				// cout << "success" << endl;
				m_matingScheme->submitScratch(pop, scratch);
				break;
			}
		}
		return true;
	}

	bool controlledBinomialSelection::mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit)
	{
		this->resetNumOffspring();
		// scrtach will have the right structure.
		this->prepareScratchPop(pop, scratch);

		size_t pldy = pop.ploidy(), nLoci=m_loci.size();
		size_t i, p;

		DBG_DO(DBG_MATING, m_famSize.clear());

		DBG_ASSERT( pop.numSubPop() == scratch.numSubPop(), SystemError,
			"Number of subpopulation can not be changed.");

		// expected frequency at each locus
		vectorf expFreq;
		PyCallFunc( m_freqFunc, "(i)", pop.gen(), expFreq, PyObj_As_Array);

		// determine expected number of alleles of each allele
		// at each subpopulation.
		UINT numSP = pop.numSubPop();

		vectoru expAlleles( nLoci * numSP );
		if( numSP > 1 && expFreq.size() == nLoci )
			// exp frequencies in subpopulation is not specified.
		{
			vectorlu totalAlleles( nLoci );
			for(i=0; i< nLoci; ++i)
			{
				int locus = m_loci[i];
				Allele allele = m_alleles[i];

				if( expFreq[i] < 0. )
					expFreq[i] = 0.;
				if( expFreq[i] > 1. )
					expFreq[i] = 1.;
				totalAlleles[i] = static_cast<ULONG>(expFreq[i]*pop.popSize()*pop.ploidy());
				// make sure one allele (our seed :-) exists
				if( expFreq[i]>0. && totalAlleles[i] == 0)
					totalAlleles[i] = 1;

				// determine the number alleles at each subpopulation.
				vectorf curFreq( numSP );
				ULONG numOfAlleles=0;
				for( size_t sp=0; sp < numSP; ++sp)
				{
					ULONG n=0;
					if(pop.shallowCopied())
					{
						for(population::IndIterator it=pop.indBegin(sp); it < pop.indEnd(sp); ++it)
							for(p=0; p<pldy; ++p)
								if( it->allele(locus, p) == allele )
									n++;
					}
					else
					{
						for(GappedAlleleIterator a=pop.alleleBegin(locus, sp),
							aEnd=pop.alleleEnd(locus, sp); a != aEnd; ++a)
						{
							if( AlleleUnsigned(*a) == allele )
								n++;
						}
					}
					numOfAlleles += n;
					curFreq[sp] = double(n)/(pop.subPopSize(sp)*pldy);
				}

				DBG_DO(DBG_MATING, cout << "Cur freq " << curFreq << endl);

				// if there is no alleles
				if( numOfAlleles == 0 && expFreq[i] > 0.)
					throw ValueError("No disease allele exists, but exp allele frequency is greater than 0.\n"
						" Generation " + toStr(pop.gen()) );

				/// calculate exp number of affected offspring in the next generation.
				///
				/// step 1: totalsize*expFreq is the total number of disease alleles
				/// step 2: assign these alleles to each subpopulation according to a multi-nomial
				/// distribution with p_i beging allele frequency at each subpopulation.
				// assign these numbers to each subpopulation
				rng().randMultinomial(static_cast<unsigned int>(pop.popSize()*expFreq[i]*pldy),
					curFreq, expAlleles.begin()+numSP*i);
			}
		}
		else if( expFreq.size() == numSP*nLoci )  // simpler case, one subpopulation, or with gieven allele frequency
		{
			for(i=0; i< nLoci; ++i)
			{
				for( size_t sp=0; sp < numSP; ++sp)
				{
#ifndef OPTIMIZED
					int locus = m_loci[i];
					Allele allele = m_alleles[i];
					ULONG n=0;
					// go through all alleles
					if(pop.shallowCopied())
					{
						for(population::IndIterator it=pop.indBegin(sp); it < pop.indEnd(sp); ++it)
							for(p=0; p<pldy; ++p)
								if( it->allele(locus, p) == allele )
									n++;
					}
					else
					{
						for(GappedAlleleIterator a=pop.alleleBegin(locus, sp),
							aEnd=pop.alleleEnd(locus, sp); a != aEnd; ++a)
						{
							if( AlleleUnsigned(*a) == allele )
								n++;
						}
					}

					// if there is no alleles
					if( n == 0 && expFreq[numSP*i+sp] > 0.)
						throw ValueError("No disease allele exists, but exp allele frequency is greater than 0.\n"
							" Generation " + toStr(pop.gen()) );
#endif
					expAlleles[numSP*i+sp] = static_cast<UINT>(pop.subPopSize(sp)*pldy*expFreq[numSP*i+sp]);
					if( expFreq[numSP*i+sp] > 0. && expAlleles[numSP*i+sp] == 0)
						expAlleles[numSP*i+sp] = 1;
				}
			}									  // each locus
		}
		else
		{
			throw ValueError("Returned expected frequency has wrong length");
		}

		DBG_DO(DBG_MATING, cout << "expected alleles " << expAlleles << endl);

		//
		GappedInfoIterator fitness;
        bool selectionOn = pop.hasVar("selection") and pop.getVarAsBool("selection");
        if (selectionOn)
        {
            UINT fit_id = pop.infoIdx("fitness");
	    	fitness = pop.infoBegin(fit_id);
        }
		/// determine if mate() will generate offspring genotype
		bool formOffGeno = this->formOffGenotype(ops);

		// for each subpopulation
		for(UINT sp=0; sp < pop.numSubPop(); ++sp)
		{
			UINT spSize = pop.subPopSize(sp);

			// total allowed disease alleles.
			vectoru totAllele( nLoci );
			// currently available disease allele.
			vectoru curAllele( nLoci, 0 );
			for(i=0; i< nLoci; ++i)
			{
				totAllele[i] = expAlleles[sp+numSP*i];
				if( totAllele[i] > spSize*pldy )
				{
					cout << "Warning: number of planned affected alleles exceed population size.";
					totAllele[i] = spSize*pldy;
				}
			}

			if( spSize == 0 )
				continue;

			// if selection is on
			if( selectionOn )
				m_sampler.set( vectorf(fitness + pop.subPopBegin(sp),
					fitness+ pop.subPopEnd(sp) ) );

			ULONG spInd = 0;
			ULONG spIndEnd = scratch.subPopSize(sp);

			// ploidy
			vectori na(nLoci, 0);
			bool allDone = true;
			for(i=0; i<nLoci; ++i)
			{
				if( totAllele[i] > 0 )
				{
					allDone = false;
					break;
				}
			}
			bool hasAff;
			int noAAcount = 0;
			// start!
			while( spInd < spIndEnd )
			{
				individual * parent;
				// choose a parent
				if( selectionOn)
					parent = &pop.ind( m_sampler.get(), sp);
				else
					parent = &pop.ind( rng().randInt(spSize), sp);

				for(i=0; i<nLoci; ++i)
					na[i] = 0;
				hasAff = false;

				// generate m_numOffspring offspring per mating
				UINT numOS=0, numOSEnd;
				for(numOS=0, numOSEnd = this->numOffspring(pop.gen() );
					numOS < numOSEnd && spInd+numOS < spIndEnd;  numOS++)
				{
					population::IndIterator it = scratch.indBegin(sp) + spInd + numOS;

					if(  formOffGeno )
						/// use deep copy!!!!!!!
						it->copyFrom(*parent);

					// apply during mating operators
					for( vector<Operator *>::iterator iop = ops.begin(),
						iopEnd = ops.end(); iop != iopEnd;  ++iop)
					{
						try
						{
							// During mating operator might reject this offspring.
							if(!(*iop)->applyDuringMating(pop, it, parent, NULL))
							{
								numOS --;
								break;
							}
						}
						catch(...)
						{
							cout << "DuringMating operator " << (*iop)->__repr__() << " throws an exception." << endl << endl;
							throw;
						}
					}							  // all during-mating operators

					// success count na, nu
					for(i=0; i<nLoci; ++i)
					{
						for(p=0; p<pldy; ++p)
						{
							if( it->allele(m_loci[i], p) == m_alleles[i] )
							{
								na[i]++;
								hasAff = true;
							}
						}
					}
				}								  // family generated.

				// now check if this family is usable.
				// later on, we can add more complicated acceptance scheme
				bool accept = false;
				// all disease alleles have been satidfied.
				// only accept unaffected families
				if( allDone )
				{
					if(!hasAff)
					{
						// has AA, so no need to compromise
						noAAcount = -1;
						accept = true;
					}
					// tried 100 times, no AA is found.
					else if( noAAcount >= 100 )
					{
						if(noAAcount == 100)
							cout << "Warning: there might not be any AA, accept unqualified individuals." << endl;
						noAAcount++;
						accept = true;
					}
					else if( noAAcount >= 0 )
						noAAcount ++;
				}
				// otherwise, accept any family that can help.
				else
				{
					for(i=0; i<nLoci; ++i)
					{
						// accept the whole family, if we need this allele
						if( curAllele[i] < totAllele[i] && na[i] > 0 )
							accept = true;
					}
				}
				// reject this family
				if( ! accept)
					continue;

				// accpet this
				if( ! allDone)
				{
					allDone = true;
					for(i=0; i<nLoci; ++i)
					{
						curAllele[i] += na[i];
						if( curAllele[i] < totAllele[i] )
							allDone = false;
					}
				}

				// accept
				DBG_DO(DBG_DEVEL, cout << curAllele << " " << spInd <<  " out of " << spIndEnd << endl);

				// success
				spInd += numOS;
			}
			if( ! allDone )
				cout << "Can not obtain enough disease alleles at generation " << pop.gen() << endl;
		}										  // all subpopulation.

		if(submit)
			submitScratch(pop, scratch);
		return true;
	}


	bool controlledRandomMating::mate( population& pop, population& scratch, vector<Operator *>& ops, bool submit)
	{
		bool hasSexChrom = pop.sexChrom();

		this->resetNumOffspring();
		// scrtach will have the right structure.
		this->prepareScratchPop(pop, scratch);

		DBG_DO(DBG_MATING, m_famSize.clear());

		DBG_ASSERT( pop.numSubPop() == scratch.numSubPop(), SystemError,
			"Number of subpopulation can not be changed.");

		size_t pldy = pop.ploidy(), nLoci=m_loci.size();
		size_t i, p;

		DBG_DO(DBG_MATING, m_famSize.clear());

		DBG_ASSERT( pop.numSubPop() == scratch.numSubPop(), SystemError,
			"Number of subpopulation can not be changed.");

		// expected frequency at each locus
		vectorf expFreq;
		PyCallFunc( m_freqFunc, "(i)", pop.gen(), expFreq, PyObj_As_Array);

		DBG_DO(DBG_MATING, cout << "expected freq " << expFreq << endl);
		// determine expected number of alleles of each allele
		// at each subpopulation.
		UINT numSP = pop.numSubPop();

		vectoru expAlleles( nLoci * numSP );
		if( numSP > 1 && expFreq.size() == nLoci )
			// exp frequencies in subpopulation is not specified.
		{
			vectorlu totalAlleles( nLoci );
			for(i=0; i< nLoci; ++i)
			{
				int locus = m_loci[i];
				Allele allele = m_alleles[i];

				if( expFreq[i] < 0. )
					expFreq[i] = 0.;
				if( expFreq[i] > 1. )
					expFreq[i] = 1.;
				totalAlleles[i] = static_cast<ULONG>(expFreq[i]*pop.popSize()*pop.ploidy());
				// make sure one allele (our seed :-) exists
				if( expFreq[i]>0. && totalAlleles[i] == 0)
					totalAlleles[i] = 1;

				// determine the number alleles at each subpopulation.
				vectorf curFreq( numSP );
				ULONG numOfAlleles=0;
				//
				for( size_t sp=0; sp < numSP; ++sp)
				{
					ULONG n=0;
					// go through all alleles
					if(pop.shallowCopied())
					{
						for(population::IndIterator it=pop.indBegin(sp); it < pop.indEnd(sp); ++it)
							for(p=0; p<pldy; ++p)
								if( it->allele(locus, p) == allele )
									n++;
					}
					else
					{
						for(GappedAlleleIterator a=pop.alleleBegin(locus, sp),
							aEnd=pop.alleleEnd(locus, sp); a != aEnd; ++a)
						{
							if( AlleleUnsigned(*a) == allele )
								n++;
						}
					}
					numOfAlleles += n;

					curFreq[sp] = double(n)/(pop.subPopSize(sp)*pldy);
				}

				// if there is no alleles
				if( numOfAlleles == 0 && expFreq[i] > 0.)
					throw ValueError("No disease allele exists, but exp allele frequency is greater than 0.\n"
						" Generation " + toStr(pop.gen()) );

				/// calculate exp number of affected offspring in the next generation.
				///
				/// step 1: totalsize*expFreq is the total number of disease alleles
				/// step 2: assign these alleles to each subpopulation according to a multi-nomial
				/// distribution with p_i beging allele frequency at each subpopulation.
				// assign these numbers to each subpopulation
				rng().randMultinomial(static_cast<unsigned int>(scratch.popSize()*expFreq[i]*pldy),
					curFreq, expAlleles.begin()+numSP*i);

				DBG_DO(DBG_MATING, cout << "DSL " << i << " Cur freq: " << curFreq << " New num "
					<< vectori(expAlleles.begin()+numSP*i, expAlleles.begin()+numSP*(i+1)) << endl);
			}
		}
		else if( expFreq.size() == numSP*nLoci )  // simpler case, one subpopulation, or with gieven allele frequency
		{
			for(i=0; i< nLoci; ++i)
			{
				for(size_t sp=0; sp < numSP; ++sp)
				{
#ifndef OPTIMIZED
					int locus = m_loci[i];
					Allele allele = m_alleles[i];
					ULONG n=0;
					// go through all alleles
					if(pop.shallowCopied())
					{
						for(population::IndIterator it=pop.indBegin(sp); it < pop.indEnd(sp); ++it)
							for(p=0; p<pldy; ++p)
								if( it->allele(locus, p) == allele )
									n++;
					}
					else
					{
						for(GappedAlleleIterator a=pop.alleleBegin(locus, sp),
							aEnd=pop.alleleEnd(locus, sp); a != aEnd; ++a)
						{
							if( AlleleUnsigned(*a) == allele )
								n++;
						}
					}

					// if there is no alleles
					if( n == 0 && expFreq[numSP*i+sp] > 0.)
						throw ValueError("No disease allele exists, but exp allele frequency is greater than 0.\n"
							" Generation " + toStr(pop.gen()) );
#endif
					expAlleles[numSP*i+sp] = static_cast<UINT>(pop.subPopSize(sp)*pldy*expFreq[numSP*i+sp] );
					if( expFreq[numSP*i+sp] > 0. && expAlleles[numSP*i+sp] == 0)
						expAlleles[numSP*i+sp] = 1;
				}
			}									  // each locus
		}
		else
		{
			throw ValueError("Returned expected frequency has wrong length");
		}

		DBG_DO(DBG_MATING, cout << "expected alleles " << expAlleles << endl);

		// empty fitness means no selection
		GappedInfoIterator fitness;
        bool selectionOn = pop.hasVar("selection") and pop.getVarAsBool("selection");
        if (selectionOn)
        {
            UINT fit_id = pop.infoIdx("fitness");
	    	fitness = pop.infoBegin(fit_id);
        }

		/// determine if any during-mating operator will generate offspring genotype
		bool formOffGeno = this->formOffGenotype(ops);

		UINT numMale, numFemale;

		/// random mating happens within each subpopulation
		for(UINT sp=0; sp < pop.numSubPop(); ++sp)
		{
			ULONG spSize = pop.subPopSize(sp);

			// total allowed disease alleles.
			vectoru totAllele( nLoci );
			// currently available disease allele.
			vectoru curAllele( nLoci, 0 );
			for(i=0; i< nLoci; ++i)
			{
				totAllele[i] = expAlleles[sp+numSP*i];
				if( totAllele[i] > spSize*pldy )
				{
					cout << "Warning: number of planned affected alleles exceed population size.";
					totAllele[i] = spSize*pldy;
				}
			}

			if( spSize == 0 ) continue;

			numMale = 0;
			for( population::IndIterator it=pop.indBegin(sp), itEnd = pop.indEnd(sp); it < itEnd;  ++it)
				if(it->sex() == Male)
					numMale ++;

			// to gain some performance, allocate memory at first.
			m_maleIndex.resize(numMale);
			m_femaleIndex.resize(spSize-numMale);

			numMale = 0;
			numFemale = 0;

			for( ULONG it=pop.subPopBegin(sp), itEnd =pop.subPopEnd(sp); it < itEnd;  it++)
			{
				if( pop.ind(it).sex() == Male)
					m_maleIndex[numMale++] = it;
				else
					m_femaleIndex[numFemale++] = it;
			}

			/// now, all individuals of needToFind sex is collected
			if( (numMale == 0 || numFemale ==0 ) && !m_contWhenUniSex )
			{
				throw ValueError("Subpopulation becomes uni-sex. Can not continue. \n"
					"You can use ignoreParentsSex (do not check parents' sex) or \ncontWhenUnixSex "
					"(same sex mating if have to) options to get around this problem.");
			}

			DBG_ASSERT( numFemale + numMale == spSize, SystemError,
				"Wrong number of male/female.");

			/// if selection is on
			if( selectionOn)
			{
				m_maleFitness.resize(numMale);
				m_femaleFitness.resize(numFemale);

				size_t ind;

				for( ind = 0; ind < numMale; ++ind)
					m_maleFitness[ind] = fitness[ m_maleIndex[ind] ];
				for( ind = 0; ind < numFemale; ++ind)
					m_femaleFitness[ind] = fitness[ m_femaleIndex[ind] ];

				m_malesampler.set(m_maleFitness);
				m_femalesampler.set(m_femaleFitness);
			}

			// generate scratch.subPopSize(sp) individuals.
			ULONG spInd = 0;
			ULONG spIndEnd = scratch.subPopSize(sp);
			// ploidy
			vectori na(nLoci, 0);
			bool allDone = true;
			for(i=0; i<nLoci; ++i)
			{
				if( totAllele[i] > 0 )
				{
					allDone = false;
					break;
				}
			}
			bool hasAff;
			// start!
			while( spInd < spIndEnd )
			{
				// randomly choose parents
				individual * dad, *mom;
				RNG& rnd = rng();

				if( selectionOn )			  // with selection
				{
					// using weidhted sampler.
					if( numMale != 0 )
						dad = &pop.ind( m_maleIndex[ m_malesampler.get() ] );
					else
						dad = &pop.ind( m_femaleIndex[ m_femalesampler.get() ] );

					if( numFemale != 0 )
						mom = &pop.ind( m_femaleIndex[ m_femalesampler.get() ] );
					else
						mom = &pop.ind( m_maleIndex[ m_malesampler.get() ]);
				}
				else
				{
					// using random sample.
					if( numMale != 0 )
						dad = &pop.ind( m_maleIndex[ rnd.randInt(numMale) ]);
					else
						dad = &pop.ind( m_femaleIndex[ rnd.randInt(numFemale) ]);

					if( numFemale != 0 )
						mom = &pop.ind( m_femaleIndex[ rnd.randInt(numFemale) ]);
					else
						mom = &pop.ind( m_maleIndex[ rnd.randInt(numMale) ]);
				}

				for(i=0; i<nLoci; ++i)
					na[i] = 0;
				hasAff = false;
				int noAAcount = 0;
				// generate m_numOffspring offspring per mating
				UINT numOS=0, numOSEnd;
				for(numOS=0, numOSEnd = this->numOffspring(pop.gen());
					numOS < numOSEnd && spInd+numOS < spIndEnd;  numOS++)
				{
					population::IndIterator it = scratch.indBegin(sp) + spInd + numOS;

					//
					// assign sex randomly
					if( ! hasSexChrom)
					{
						int offSex = rnd.randInt(2);
						it->setSex( offSex==0?Male:Female);
					}

					if( formOffGeno )			  // use the default no recombination random mating.
					{
												  // initialize to avoid compiler complains
						int dadPloidy=0, momPloidy=1;

						for(UINT ch=0, chEnd = dad->numChrom(); ch < chEnd;  ++ch)
						{
							dadPloidy = rnd.randInt(2);
							momPloidy = rnd.randInt(2);

							DBG_ASSERT((dadPloidy==0 || dadPloidy==1) &&
								( momPloidy==0 || momPloidy==1), ValueError,
								"Ploidy must be 0 or 1");

							copy(dad->genoBegin(dadPloidy, ch),
								dad->genoEnd(dadPloidy,ch) ,
								it->genoBegin(0,ch));
							copy(mom->genoBegin(momPloidy,ch),
								mom->genoEnd(momPloidy,ch) ,
								it->genoBegin(1,ch));
						}

						// last chromosome (sex chromosomes)
						if( hasSexChrom )
						{
							if( dadPloidy == 1)	  // Y of XY
								it->setSex(Male);
							else
								it->setSex(Female);
						}
					}

					/// apply all during mating operators
					for( vector<Operator *>::iterator iop = ops.begin(), iopEnd = ops.end(); iop != iopEnd;  ++iop)
					{
						try
						{
							// During mating operator might reject this offspring.
							if(!(*iop)->applyDuringMating(pop, it, dad, mom))
							{
								numOS--;
								break;
							}
						}
						catch(...)
						{
							cout << "DuringMating operator " << (*iop)->__repr__() << " throws an exception." << endl << endl;
							throw;
						}
					}
					// success count na, nu
					for(i=0; i<nLoci; ++i)
					{
						for(p=0; p<pldy; ++p)
						{
							if( it->allele(m_loci[i], p) == m_alleles[i] )
							{
								na[i]++;
								hasAff = true;
							}
						}
					}

				}								  // family generated.

				// now check if this family is usable.
				// later on, we can add more complicated acceptance scheme
				bool accept = false;
				// all disease alleles have been satidfied.
				// only accept unaffected families
				// otherwise, accept any family that can help.
				if( allDone )
				{
					if(!hasAff)
					{
						// has AA, so no need to compromise
						noAAcount = -1;
						accept = true;
					}
					// tried 100 times, no AA is found.
					else if( noAAcount >= 100 )
					{
						if(noAAcount == 100)
							cout << "Warning: there might not be any AA, accept unqualified individuals." << endl;
						noAAcount++;
						accept = true;
					}
					else if( noAAcount >= 0 )
						noAAcount ++;
				}
				else
				{
					// acceptance scheme
					if( m_acceptScheme == CM_AcceptOneDSL)
					{
						for(i=0; i<nLoci; ++i)
						{
							// accept the whole family, if we need this allele
							if( curAllele[i] < totAllele[i] && na[i] > 0 )
								accept = true;
						}
					}
					else if( m_acceptScheme == CM_AcceptAllDSL)
					{
						if(hasAff)
							accept = true;
					}
					else
						throw ValueError("UNrecognized acceptance scheme");
				}
				// reject this family
				if( ! accept)
					continue;

				// accpet this family
				if( ! allDone)
				{
					allDone = true;
					for(i=0; i<nLoci; ++i)
					{
						curAllele[i] += na[i];
						if( curAllele[i] < totAllele[i] )
							allDone = false;
					}
				}

				// accept
				DBG_DO(DBG_DEVEL, cout << curAllele << " " << spInd <<  " out of " << spIndEnd << endl);

				// success
				spInd += numOS;

				DBG_WARNING( spInd == spIndEnd && (! allDone),
					cout << "Offspring generated finished with allele freq unsatisfied." << endl);
			}
		}										  // each subPop

		if(submit)
			submitScratch(pop, scratch);
		return true;
	}

}
