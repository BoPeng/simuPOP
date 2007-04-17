/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
 *                                                                         *
 *   $LastChangedDate: 2006-02-21 15:27:25 -0600 (Tue, 21 Feb 2006)        *
 *   $Rev: 191$                                                            *
 *                                                                         *
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

#include "misc.h"

namespace simuPOP
{
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
		DBG_FAILIF( curGen == 0 && (ValidPyObject(NtFunc) || ValidPyObject(fitnessFunc)),
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
		{
			// default to [1,1,1]
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

				// for extremely small a, assume a=0,
				if( fabs(a) < 1e-8 )
				{
					y1 = -c/b;
					// y1 should be valid
					y2 = 1000.;
				}
				else
				{
					y1 = (-b+sqrt(b2_4ac))/(2*a);
					y2 = (-b-sqrt(b2_4ac))/(2*a);
				}

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
			{
				// when the allele get fixed, restart
				if(allowFixation)
				{
					break;
				}
				else
				{
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

	//
	// for example
	//      BB  Bb   bb
	//  AA  0   1    2
	//  Aa  3   4    5
	//  aa  6   7    8
	//
	// allgeno = [1, 2] = Aa, bb
	// first loop
	//      index = 0*3+1
	// second loop
	//      index = 1*3 + 2 = 5
	/// CPPONLY
	double fitOfGeno(unsigned loc, const vectori & allgeno, const vectorf & fitness, const vectorf::const_iterator & freq)
	{
		int index = 0;
		double fq = 1;
		for(size_t i=0; i<allgeno.size(); ++i)
		{
			if(i != loc)
			{
				if(allgeno[i] == 0)
					fq *= (1-freq[i])*(1-freq[i]);
				else if(allgeno[i] == 1)
					fq *= 2*(1-freq[i])*freq[i];
				else
					fq *= freq[i]*freq[i];
			}
			index = index*3 + allgeno[i];
			if(fq == 0.)
				return 0.;
		}
		return fitness[index]*fq;
	}

	/// CPPONLY
	// get individual fitness, accounting interaction
	void interFitness(unsigned nLoci, const vectorf & fitness, const vectorf::const_iterator & freq, vectorf & sAll)
	{
		sAll.resize(3*nLoci);
		// each locus
		for(size_t loc=0; loc<nLoci; ++loc)
		{
			// each genotype AA, Aa and aa,
			// geno is actually the number of disease allele
			for(size_t geno=0; geno<3; ++geno)
			{
				// iterate through OTHER DSL
				vectori allgeno(nLoci, 0);
				// set myself.
				allgeno[loc] = geno;
				// iterator through others
				double f = 0;
				for(size_t it = 0; it < std::pow3(nLoci-1); it++)
				{
					// assign allgeno
					size_t num = it;
					for(size_t l=0; l<nLoci; ++l)
					{
						if(l==loc)
							continue;
						allgeno[l] = num % 3;
						num /= 3;
					}
					f += fitOfGeno(loc, allgeno, fitness, freq);
					DBG_DO(DBG_DEVEL, cout << allgeno << " " << fitOfGeno(loc, allgeno, fitness, freq) << " , " );
				}
				// sum over other genotype
				DBG_DO(DBG_DEVEL, cout << loc << " " << geno << " " <<  f << endl);
				sAll[loc*3+geno] = f;
			}
			// convert to the form 1, s1, s2
			sAll[3*loc+1] = sAll[3*loc+1] / sAll[3*loc] - 1.;
			sAll[3*loc+2] = sAll[3*loc+2] / sAll[3*loc] - 1.;
			sAll[3*loc] = 0.;
		}
		DBG_DO(DBG_MATING, cout << "fitness " << fitness << " freq " << freq[0] << ", " << " sall " << sAll << endl);
	}

	// the python version of interaction Fitness, for convenience purpose.
	/// CPPONLY
	vectorf MarginalFitness(unsigned nLoci, const vectorf & fitness, const vectorf & freq)
	{
		vectorf sAll(3, 1);
		interFitness(nLoci, fitness, freq.begin(), sAll);
		return sAll;
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
		DBG_FAILIF( curGen == 0 && (NtFunc != NULL || ValidPyObject(fitnessFunc)),
			ValueError, "curGen should be > 0 if NtFunc or fitnessFunc is defined.");
		DBG_FAILIF( curGen > 0 && curGen < maxMutAge, ValueError,
			"curGen should be >= maxMutAge");

		matrix result(nLoci);

		// in the cases of independent and constant selection pressure
		// easy case.
		if( InvalidPyObject(fitnessFunc) && fitness.size() == nLoci*3)
		{
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
		bool interaction = false;
		if(ValidPyObject(fitnessFunc))
		{
			if( !PyCallable_Check(fitnessFunc) )
				throw ValueError("sFunc is not a valid Python function.");
			else
				// increase the ref, just to be safe
				Py_INCREF(fitnessFunc);
		}
		else if( fitness.empty() )
		{
			// default to [1,1,1]
			sAll.resize(nLoci*3, 0.);
		}
		// independent case
		else if( fitness.size() == 3*nLoci)
		{
			for(i=0; i<nLoci; ++i)
			{
				// convert to the form 1, s1, s2
				sAll[3*i+1] = fitness[3*i+1] / fitness[3*i] - 1.;
				sAll[3*i+2] = fitness[3*i+2] / fitness[3*i] - 1.;
				sAll[3*i] = 0.;
			}
		}
		// interaction case
		else if (fitness.size() == std::pow3(nLoci))
		{
			interaction = true;
		}
		else
			throw ValueError("Wrong size of fitness vector");

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
			if(ValidPyObject(fitnessFunc))
			{
				// compile allele frequency... and pass
				vectorf sAllTmp;
				PyObject* freqObj = Double_Vec_As_NumArray( xt.begin()+nLoci*idx, xt.begin()+nLoci*(idx+1) );
				PyCallFunc2(fitnessFunc, "(iO)", curGen-idx-1, freqObj, sAllTmp, PyObj_As_Array);

				if(sAllTmp.size() == 3*nLoci)
				{
					for(i=0; i<nLoci; ++i)
					{
						// convert to the form 1, s1, s2
						sAll[3*i+1] = sAllTmp[3*i+1] / sAllTmp[3*i] - 1.;
						sAll[3*i+2] = sAllTmp[3*i+2] / sAllTmp[3*i] - 1.;
						sAll[3*i] = 0.;
					}
				}
				// interaction case.
				else if(sAllTmp.size() == std::pow3(nLoci))
				{
					// xt is the current allele frequency
					interFitness(nLoci, sAllTmp, xt.begin()+(idx*nLoci), sAll);
				}
				else
					throw ValueError("Wrong size of fitness vector: " + toStr(sAll.size()));
			}
			else if (interaction)
			{
				// from fitness vector, get sAll using allele frequency
				interFitness(nLoci, fitness, xt.begin()+(idx*nLoci), sAll);
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

					if( fabs(a) < 1e-8 )
					{
						y1 = -c/b;
						// y1 should be valid
						y2 = 1000.;
					}
					else
					{
						y1 = (-b+sqrt(b2_4ac))/(2*a);
						y2 = (-b-sqrt(b2_4ac))/(2*a);
					}

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
		if(ValidPyObject(fitnessFunc))
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

#ifndef OPTIMIZED
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
#endif
}
