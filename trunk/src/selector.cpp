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

#include "simupop_cfg.h"
#include "utility.h"

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
      if( ! m_phase && a > b )                    // ab=ba
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
    for(vectoru::iterator loc=m_loci.begin(); loc!=m_loci.end(); ++loc)
    {
      /// get genotype of ind
      Allele a = ind->allele(*loc, 0);
      Allele b = ind->allele(*loc, 1);

      int numWildtype=0;

      // count number of wildtype
      if( find(m_wildtype.begin(), m_wildtype.end(), a) != m_wildtype.end() )
        numWildtype ++;

      if( find(m_wildtype.begin(), m_wildtype.end(), b) != m_wildtype.end() )
        numWildtype ++;

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

    return 0.0;
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
    catch(...)                                    // if there is no sample
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
    pop.setIndInfoWithSubPopID();

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
        pop.ind( pick[i] ).setInfo( -1 );
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
          pop.ind( pick[i], sp ).setInfo( -1 );
      }                                           // each subpop
    }                                             // else
    return pop.newPopByIndInfo(false);
  }

  bool caseControlSample::prepareSample(population& pop )
  {
    if( ! m_spSample )                            // sample from the whole population.
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
    if( ! m_spSample)                             // draw sample from the whole population
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
        pop.ind( m_caseIdx[0][i] ).setInfo( 0 );
      }
      // remove others
      for( i= nCase; i < numAffected; ++i)
        pop.ind( m_caseIdx[0][i] ).setInfo( -1 );

      // keep first m_size individuals of shuffled indices
      for( i=0; i < nControl; ++i)
      {
        nControlInSP[pop.subPopIndPair(m_controlIdx[0][i]).first]++;
        pop.ind( m_controlIdx[0][i] ).setInfo( 1 );
      }
      // remove others
      for( i= nControl; i < numUnaffected; ++i)
        pop.ind( m_controlIdx[0][i] ).setInfo( -1 );

      DBG_DO(DBG_SELECTOR, cout << "Getting sample population" << endl);
      population& sample = pop.newPopByIndInfo(false);

      sample.setIntVectorVar("nCases", nCaseInSP);
      sample.setIntVectorVar("nControls", nControlInSP);
      // determine exactly how many cases and controls in the final sample
      return sample;
    }
    else                                          // sample from each subpop
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
          pop.ind( m_caseIdx[sp][i],sp ).setInfo( 0 );
        // remove others
        for(int i= nCase, iEnd= m_caseIdx[sp].size(); i<iEnd; ++i)
          pop.ind( m_caseIdx[sp][i],sp ).setInfo( -1 );

        // keep first m_size individuals of shuffled indices
        nControl = std::min(m_numControls[sp], static_cast<int>(m_controlIdx[sp].size()));
        for(int i=0; i < nControl; ++i)
          pop.ind( m_controlIdx[sp][i],sp ).setInfo( 1 );
        // remove others
        for(int i= nControl, iEnd = m_controlIdx[sp].size(); i<iEnd; ++i)
          pop.ind( m_controlIdx[sp][i],sp ).setInfo( -1 );
      }
      // newPop .... but ignore ancestral populations
      population& sample = pop.newPopByIndInfo(false);
      sample.setIntVectorVar("nCases", m_numCases);
      sample.setIntVectorVar("nControls", m_numControls);
      return sample;
    }
  }

  bool affectedSibpairSample::prepareSample(population& pop)
  {
    // get tag info for each subpop
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

      // get tag and affected status info
      std::vector<TagType> tags(spSize);
      std::vector<bool> affected(spSize);
      for(size_t i=0; i < spSize; ++i)
      {
        tags[i] = pop.ind(i, sp).tag();
        affected[i] = pop.ind(i, sp).affected();
      }

      // find affected sibpairs, difficult...
      TagType tag;
      for(size_t i=0; i< spSize; ++i)
      {
        if( affected[i] == m_affectedness )
        {
          tag = tags[i];
          // ignore father = monther which is used when no opposite
          // sex exist when mating.
          if(tag.first == tag.second)
            continue;
          // find an individual with the same parents
          for(size_t j = i+1; j < spSize; ++j)
          {
            // ok, find it
            if( affected[j] == m_affectedness && tags[j] == tag )
            {
              // use absolute index, so spBegin+...
              sibpairs.push_back(
                std::pair<ULONG,ULONG>(spBegin+i, spBegin+j) );
              // mark location i, j as occupied.
              affected[i] = ! m_affectedness;
              affected[j] = ! m_affectedness;
              // excluding indiviudus having one of the same parent.
              for(size_t k = i+1; k < spSize; ++k)
              {
                if( tags[k].first == tag.first || tags[k].second == tag.second)
                  affected[k] = ! m_affectedness;
              }
              break;
            }
          }                                       // j
        }                                         // affected [i]
      }                                           // for all i

      DBG_DO(DBG_SELECTOR, cout << "Sibpairs at SP " << sp
        << " is " << sibpairs.size() << endl);

      // now, we know the number of affected sibpairs in subpop i
      pop.setIntVar( subPopVar_String(sp, "numAffectedSibpairs"),
        sibpairs.size());

      nSibs += sibpairs.size();

      if( m_size.size() > 1 && m_size[sp] > sibpairs.size())
        cout << "Warning: Not enough sibpairs (" << sibpairs.size()
          << ") to be sampled at subpop " << sp << endl;
    }                                             // each subpop

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
    vectori chosenOff, chosenPar;
    vectori nSibpairSP(pop.numSubPop());

    if( m_size.size() <= 1 )                      // draw from the whole population
    {
      DBG_DO(DBG_SELECTOR, cout << "Draw from the whole population." << endl);

      /// sibs when m_size.size() <= 1
      vector< std::pair<ULONG, ULONG> > allSibs;

      for(UINT sp =0; sp < pop.numSubPop(); ++sp)
        allSibs.insert( allSibs.end(), m_sibpairs[sp].begin(), m_sibpairs[sp].end());

      if( !m_size.empty() && m_size[0] > allSibs.size() )
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
        ind->setInfo(-1);
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
        chosenOff.push_back( (int)(allSibs[ idx[i] ].first) );
        chosenOff.push_back( (int)(allSibs[ idx[i] ].second) );
        // DBG_ASSERT( pop.ind( allSibs[ idx[i] ].first ).info()==-1,
        //  SystemError, "Duplicate selection");
        //DBG_ASSERT( pop.ind( allSibs[ idx[i] ].second).info()==-1,
        //  SystemError, "Duplicate selection");
        pop.ind( allSibs[ idx[i] ].first ).setInfo(i);
        pop.ind( allSibs[ idx[i] ].second ).setInfo(i);
        TagType tag = pop.ind( allSibs[ idx[i]].first ).tag();
        chosenPar.push_back( (int)(tag.first) );
        chosenPar.push_back( (int)(tag.second) );
        parents.push_back(
          std::pair<ULONG, ULONG>((ULONG)(tag.first), (ULONG)(tag.second))
          );
      }

      // keep the order and count number of selected ones in each subpop
      // here we say number of individuals chosen from a sp += 2
      for( size_t i=0; i<N; ++i)
        nSibpairSP[ pop.subPopIndPair( allSibs[idx[i]].first).first] += 2;
      //
      DBG_DO(DBG_SELECTOR, cout << "#ind from population " << nSibpairSP << endl);
    }
    else                                          // for each subpop
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
          ind->setInfo(-1);

        UINT N = std::min(asSize, static_cast<size_t>(m_size[sp]));

        // sample sibpairs
        random_shuffle(idx.begin(), idx.end());
        std::sort( idx.begin(), idx.begin() + N);

        // set sibpairs
        for( size_t i=0; i< N; ++i)
        {
          chosenOff.push_back( (int)(sibpairs[ idx[i]].first) );
          chosenOff.push_back( (int)(sibpairs[ idx[i]].second) );
          // sibpairs uses absolute indices, so no sp parameter is rquired.
          pop.ind( sibpairs[ idx[i]].first ).setInfo(sibID);
          pop.ind( sibpairs[ idx[i]].second ).setInfo(sibID);
          TagType tag = pop.ind(sibpairs[ idx[i]].first ).tag();
          chosenPar.push_back( (int)(tag.first) );
          chosenPar.push_back( (int)(tag.second) );
          parents.push_back(
            std::pair<ULONG, ULONG>((ULONG)(tag.first), (ULONG)(tag.second))
            );
          sibID ++;
        }
        nSibpairSP[sp] = N*2;
      }                                           // sp
    }                                             // else
    // this is offspring population with copy of ancestral pop
    // (true means keepAncestralPops)
    population & newPop = pop.newPopByIndInfo(true);

    newPop.setIntVectorVar("numOffspring", nSibpairSP);
    newPop.setIntVectorVar("chosenOffspring", chosenOff);
    newPop.setIntVectorVar("chosenParents", chosenPar);

    DBG_DO(DBG_SELECTOR, cout << "Offspring selection done." << endl);

    vectorlu tmp(pop.numSubPop());
    for(size_t i=0; i< tmp.size(); ++i)
      tmp[i] = nSibpairSP[i];

    // do not allow to change newPop size. ONLY subpop structure.
    newPop.setSubPopStru(tmp, false);
    // set parental info
    newPop.setAncestralDepth(1);
    newPop.useAncestralPop(1);

    DBG_DO(DBG_SELECTOR, cout << "Working on parent generations." << endl);

    // set everyone to -1 (remove them)
    for(population::IndIterator ind=newPop.indBegin();
      ind != newPop.indEnd(); ++ind)
    {
      ind->setInfo(-1);
    }

    DBG_DO(DBG_SELECTOR, cout << parents.size()*2 << " parents will be kept. " << endl);

    for(size_t i=0; i< parents.size(); ++i)
    {
      newPop.ind( parents[i].first).setInfo(i);
      newPop.ind( parents[i].second).setInfo(i);
    }
    DBG_DO(DBG_SELECTOR, cout << "Setting info done." << endl);

    newPop.setSubPopByIndInfo();
    // do not allow to change newPop size. ONLY subpop structure.
    DBG_DO( DBG_SELECTOR, cout << "Set structure " << tmp << endl);
    newPop.setSubPopStru(tmp, false);
    newPop.useAncestralPop(0);
    return newPop;
  }

}
