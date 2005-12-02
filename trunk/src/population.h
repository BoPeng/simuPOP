/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$                                                     *
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
#ifndef _POPULATION_H
#define _POPULATION_H

/**
\file
\brief head file of class Population
*/

#include <vector>
#include <numeric>
using std::vector;
using std::accumulate;
using std::pair;

#include <functional>
using std::equal_to;

#include "utility.h"
#include "individual.h"

#include <fstream>
using std::ifstream;
using std::ofstream;

// used to save history population
// 0 (first parental) 1, ...., n
#include <deque>
using std::deque;

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#ifndef __NO_XML_SUPPORT__
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#endif
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
using boost::serialization::make_nvp;

namespace simuPOP
{

  /** \brief a collection of individuals with subPopulation structure

  Please refer to user's Guide for details about this object.

  */
  template< class Ind>
    class Population : public GenoStruTrait
  {
    public:

      /// expose individual type (template parameter)
      typedef Ind IndType;

      /// individual itertor, used to iterate all individuals.
      typedef typename vector<IndType>::iterator IndIterator;

      /// allele iterator, used to iterate the same allele across
      /// all ploidy and individuals.
      typedef GappedAlleleIterator AlleleIterator;

      /// simple allele, used to iterate through all alleles one by one,
      /// ignoring population structure.
      typedef typename vector<Allele>::iterator SimpleIterator;

      /// iterate through all alleles of an individual.
      typedef typename IndType::GenoIterator GenoIterator;

      /// individual tag type. (template parameter for Individual )
      typedef typename IndType::TagType TagType;

      /// individual info type
      typedef typename IndType::InfoType InfoType;

      /** @name  constructors and destructor */
      //@{

      /// create a population object with given size and genotypic structure
      /**
      \param size population size. Can be ignored if subPop is specified.
         In that case, size is sum of subPop. Default to 0.
      \param ploidy number of sets of chromosomes. Default to 2 (diploid).
      \param loci an array of numbers of loci on each chromosome. If
         not specified, assume a single locus on one chromosome. Number
         of chromosomes is determined by the size of this array.
      \param lociPos an array of loci distance for each locus. You can
         also use a nested array to specify loci distance for each chromosome.
         ( [1,2,3,4,5] or [[1,2],[3,4,5]] are both allowed for loci=[2,3])
         The default values are 1, 2, etc. on each chromosome.
      \param subPop an array of subPopulation sizes. Default value is [size]
      which means a single subpopulation of the whole population. If both
      size and subPop are given, sum of subPop should agree with size.
      \param ancestralDepth number of ancestral populations to keep. Default to 0,
      meaning only current generation will be available. If -1 is given,
      all ancestral populations will be saved. This may exhaust your RAM
      pretty quickly though.
      \param alleleNames an array of allele names. The first element should be
      given to invalid/unknown allele. For example, for a locus with alleles
      A,C,T,G, you can specify alleleName as ('_','A','C','T','G'). Note that
      simuPOP uses 1,2,3,4 internally and these names will only be used for
      display purpose.
      \param lociNames an array or a matrix (separated by chromosome) of names for
      each loci. Default to "locX-X" where
      X-X is chromosome-loci index starting from 1. This info is rarely used.
      \param maxAllele maximum allele number. Default to the max allowed allele states
      of current library (standard or long allele version)
      \return no return value. Exception will be thrown is wrong parameters are given.
      \sa simulator, baseOperator, mating schemes
      \test popInit.log \include popInit.log
      */
      Population( ULONG size=0,
        UINT ploidy=2,
        const vectoru& loci=vectoru(),
        bool sexChrom=false,
        const vectorf& lociPos=vectorf(),
        const vectorlu& subPop=vectorlu(),
        int ancestralDepth=0,
        const vectorstr& alleleNames=vectorstr(),
        const vectorstr& lociNames=vectorstr(),
        UINT maxAllele = MaxAllele )
        :
      GenoStruTrait(),
        m_popSize(size),
        m_numSubPop(subPop.size()),
        m_subPopSize(subPop),
        m_popGenoSize(0),
        m_subPopIndex(subPop.size()+1),
        m_genotype(0),                            // resize later
        m_inds(0),                                // default constructor will be called.
        m_ancestralDepth(ancestralDepth),
        m_vars(NULL, true),                       // invalid shared variables initially
        m_ancestralPops(0),                       // no history first
        m_rep(-1),
        m_grp(-1),
        m_gen(0),
        m_curAncestralPop(0)
      {
        DBG_FAILIF(maxAllele > MaxAllele, ValueError,
          "maxAllele is bigger than maximum allowed allele state of this library (" + toStr(MaxAllele) +
          ")\nPlease use simuOpt.setOptions(longAllele=True) to use the long allele version of simuPOP.");

        DBG_DO( DBG_POPULATION, cout << "Constructor of Population is called\n");

        // if specify subPop but not m_popSize
        if( !subPop.empty() )
        {
          if( size == 0 )
            m_popSize = accumulate(subPop.begin(), subPop.end(), 0UL);
          else
            DBG_ASSERT( m_popSize == accumulate(subPop.begin(), subPop.end(), 0UL),
              ValueError, "If both size and subPop are specified, size should equal to sum(subPop)");
        }

        // get a GenoStructure with parameters. GenoStructure may be shared by some populations
        // a whole set of functions ploidy() etc in GenoStruTriat can be used after this step.
        this->setGenoStructure(ploidy, loci, sexChrom, lociPos, alleleNames, lociNames, maxAllele );

        DBG_DO( DBG_DEVEL, cout << "Individual size is " << sizeof(IndType) << '+'
          << sizeof(Allele) << '*' << genoSize() << endl
          << "Info: " << sizeof(InfoType) << ", Tag: " << sizeof(TagType)
          << ", GenoPtr: " << sizeof(Allele*) << ", Flag: " << sizeof(unsigned char)
          << ", plus genoStru" << endl );

        // size of genotypic data for the whole population
        m_popGenoSize = totNumLoci() * ploidy * m_popSize;

        try
        {
          // allocate memory here (not in function definition)
          m_inds.resize(m_popSize);

          // create genotype vector holding alleles for all individuals.
          m_genotype.resize( m_popGenoSize);

          // set subPopulation indexes, do not allow popsize change
          setSubPopStru(subPop, false);

          // set individual pointers
          // reset individual pointers
          Allele* ptr = &(m_genotype[0]);
          UINT step = genoSize();
          for(ULONG i=0; i< m_popSize; ++i, ptr+=step)
          {
            m_inds[i].setGenoPtr( ptr );
            m_inds[i].setGenoStruIdx(this->genoStruIdx());
            m_inds[i].setShallowCopied(false);
          }
          IndType::clearShallowCopiedFlag();
        }
        catch(...)
        {
          cout << "You are creating a population of size " << m_popSize << endl
            << " which requires approximately " << static_cast<double>(m_popGenoSize)/1024/1024
            + static_cast<double>(m_popSize) * sizeof( IndType ) /1024/1024 << "M RAM." << endl;
          cout << "Memory allocation fail. A population of size 1 is created." << endl;
          *this = Population(0);
          throw OutOfMemory("Memory allocation fail");
        }
        // set local variable
        setRep(-1);
        setGrp(-1);
      }

      /// CPPONLY copy constructor
      Population(const Population& rhs):
      GenoStruTrait(rhs),
        m_popSize(rhs.m_popSize),
        m_numSubPop(rhs.m_numSubPop),
        m_subPopSize(rhs.m_subPopSize),
        m_popGenoSize(rhs.m_popGenoSize),
        m_subPopIndex(rhs.m_subPopIndex),
        m_genotype(0),
        m_inds(0),
        m_ancestralDepth(rhs.m_ancestralDepth),
        m_vars(rhs.m_vars),                       // variables will be copied
        m_rep(-1),                                // rep is set to -1 for new pop (until simulator really set them
        m_grp(-1),
        m_gen(0),
        m_curAncestralPop(rhs.m_curAncestralPop)
      {
        DBG_DO(DBG_POPULATION,
          cout << "Copy constructor of Population is called\n" << endl);

        try
        {
          m_inds.resize(rhs.m_popSize);
          m_genotype.resize(rhs.m_popGenoSize);
        }
        catch(...)
        {
          cout << "You are creating a population of size " << m_popSize << endl
            << " which requires approximately " << static_cast<double>(m_popGenoSize)/1024/1024
            + static_cast<double>(m_popSize) * sizeof( IndType ) /1024/1024 << "M RAM." << endl;
          cout << "Memory allocation fail. A population of size 1 is created." << endl;
          *this = Population(0);
          throw OutOfMemory("Memory allocation fail");
        }

        // individuals will always have the correct genostructure
        // by using their copied pointer
        // population, however, need to set this pointer correctly
        //
        setGenoStruIdx(rhs.genoStruIdx());

        // copy genotype one by one so individual genoPtr will not
        // point outside of subPopulation region.
        Allele * ptr = &(m_genotype[0]);
        UINT step = this->genoSize();
        for(ULONG i=0; i< m_popSize; ++i, ptr+=step)
        {
          m_inds[i].setGenoPtr( ptr );
          m_inds[i].copyFrom( rhs.m_inds[i]);
        }
        IndType::clearShallowCopiedFlag();

        // copy ancestral populations
        try
        {
          // copy all. individual will be shallow copied
          m_ancestralPops = rhs.m_ancestralPops;
          // need to setGenoPtr
          for( size_t ap = 0; ap < m_ancestralPops.size(); ++ap)
          {
            popData& lp = m_ancestralPops[ap];
            const popData& rp = rhs.m_ancestralPops[ap];
            vector<IndType>& linds = lp.m_inds;
            const vector<IndType>& rinds = rp.m_inds;
            Allele* lg = &(lp.m_genotype[0]);
            const Allele* rg = &(rp.m_genotype[0]);
            ULONG ps = rinds.size();

            for(ULONG i=0; i<ps; ++i)
              linds[i].setGenoPtr( rinds[i].genoPtr() - rg + lg );
          }
        }
        catch(...)
        {
          cout << "Unable to copy ancestral populations. "
            << "The popolation size may be too big." << endl
            << "The population will still be usable but without any ancestral population stored." << endl;
          m_ancestralDepth = 0;
          m_ancestralPops.clear();
        }

        // set local variable
        setRep(-1);
        setGrp(-1);
      }

      ///
      Population& clone(bool keepAncestralPops=true) const
      {
        Population& p = *new Population<IndType>(*this);
        if( keepAncestralPops == false)
          p.m_ancestralPops.clear();
        return p;
      }

      /// SWAP population
      /// swap the content of two populations
      void swap(Population& rhs)
      {
        GenoStruTrait::swap(rhs);
        std::swap(m_popSize, rhs.m_popSize);
        std::swap(m_numSubPop, rhs.m_numSubPop);
        m_subPopSize.swap(rhs.m_subPopSize);
        std::swap(m_popGenoSize, rhs.m_popGenoSize);
        m_subPopIndex.swap(rhs.m_subPopIndex);
        m_genotype.swap(rhs.m_genotype);
        m_inds.swap(rhs.m_inds);
        std::swap(m_ancestralDepth, rhs.m_ancestralDepth);
        m_vars.swap(rhs.m_vars);
        m_ancestralPops.swap(rhs.m_ancestralPops);
        std::swap(m_rep, rhs.m_rep);
        std::swap(m_grp, rhs.m_grp);
        std::swap(m_gen, rhs.m_gen);
        std::swap(m_curAncestralPop, rhs.m_curAncestralPop);
      }

      /// destroy a population
      ~Population()
      {

        DBG_DO(DBG_POPULATION,
          cout << "Destructor of Population is called" << endl);

        // geno structure is shared so will not be removed.
      }

      // allow str(population) to get something better looking
      string __repr__()
      {
        return "<simuPOP::Population of size " + toStr(popSize()) + ">";
      }

      /// \brief set population/subPopulation given subpopulation sizes
      ///
      /// CPPONLY
      /// \param subPopSize an array of subpopulation sizes
      ///    the population may or may not change according to
      ///    parameter allowPopSizeChange if sum of subPopSize
      ///    does not match popSize.
      /// \param allowPopSizeChange if true, popSize can change to sum of
      ///    subPopSize.
      ///
      /// \return none
      /// \sa migration, mating
      /// \note this function will be used when setting up new population
      ///   during the creation of population or during mating (creation
      ///   of scratch population), users are not supposed to call it
      ///   directly. (that is why CPPONLY is set)
      ///
      void setSubPopStru(const vectorlu& subPopsize, bool allowPopSizeChange)
      {
        // case 1: remove all subpopulation structure
        // do not change population size
        // individuals are valid....
        if ( subPopsize.empty() )
        {
          m_numSubPop = 1;
          m_subPopSize.resize(1, m_popSize);
          m_subPopIndex.resize(2);
        }
        else                                      // may change populaiton size
        {
          m_numSubPop = subPopsize.size();
          m_subPopSize = subPopsize;
          m_subPopIndex.resize( m_numSubPop+1);

          ULONG totSize = accumulate(subPopsize.begin(), subPopsize.end(), 0UL);

          // usually, totSize == m_popSize, individuals are valid
          if( totSize != m_popSize)
          {
            DBG_DO( DBG_POPULATION, "Populaiton size changed to " + toStr(totSize) +
              " Genotype information may be lost");

#ifndef OPTIMIZED
            if( !allowPopSizeChange)
            {
              cout << "Total new size " << totSize << endl;
              cout << "Attempted subpop " << subPopsize << endl;
              cout << "Total current " << m_popSize << endl;
              cout << "Current subpop size " << this->subPopSizes() << endl;
              throw ValueError("Populaiton size is fixed (by allowPopSizeChange parameter)."
                " Subpop sizes should add up to popsize");
            }
#endif

            // change populaiton size
            // genotype and individual info will be kept
            // but pointers need to be recalibrated.
            m_popSize = totSize;
            m_popGenoSize = genoSize() * m_popSize;

            try
            {
              m_genotype.resize( m_popGenoSize );
              m_inds.resize(m_popSize);
            }
            catch(...)
            {
              cout << "You are creating a population of size " << m_popSize << endl
                << " which requires approximately " << static_cast<double>(m_popGenoSize)/1024/1024
                + static_cast<double>(m_popSize) * sizeof( IndType ) /1024/1024 << "M RAM." << endl;
              cout << "Memory allocation fail. " << endl;
              throw OutOfMemory("Memory allocation fail");
            }
            // reset individual pointers
            Allele* ptr = &m_genotype[0];
            UINT step = genoSize();
            for(ULONG i=0; i< m_popSize; ++i, ptr+=step)
            {
              m_inds[i].setGenoPtr( ptr );
              m_inds[i].setGenoStruIdx(genoStruIdx());
              m_inds[i].setShallowCopied(false);
            }
            IndType::clearShallowCopiedFlag();
          }
        }

        // build subPop index
        UINT i = 1;
        for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
          m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
      }

      ///  number of sub populations.
      /**
       \return number of subPopulations (>=1)
       */
      UINT numSubPop() const
      {
        return m_numSubPop;
      }

      /// get size of subPopulation subPop
      /** \param subPop index of subPopulation (start from 0)
          \return size of subPopulation subPop
      */
      ULONG subPopSize(UINT subPop) const
      {
        CHECKRANGESUBPOP(subPop);

        return m_subPopSize[subPop];
      }

      /// get size of all subPopulations
      /** \return an array of size of subPopulations
       */
      vectorlu subPopSizes()
      {
        return m_subPopSize;
      }
      //@}

      /** @name indices
        conversion between absoluate indices and relative indices.
        return of chromomosome/subPopulation indices.
      */
      //@{

      /// get population size
      /**
       \return total number of individuals in this population
       */
      ULONG popSize() const
      {
        return m_popSize;
      }

      ///  absolute index of individual at a subPopulation
      /**
        \param index index of individual at subPopulation subPop
        \param subPop subpopulation index
        \return absolute index of individual at subPopulation subPop
        \sa subPopIndPair
      */
      ULONG absIndIndex(ULONG index, UINT subPop) const
      {
        CHECKRANGESUBPOP(subPop);
        CHECKRANGESUBPOPMEMBER(index, subPop);

        return( m_subPopIndex[subPop] + index);
      }

      /// subPop and relative index of an individual
      /*
        \param absInd absolute index of an individual
        \return a pair of values (subPop, index)
        \sa absIndIndex
      */
      std::pair<UINT, ULONG> subPopIndPair(ULONG ind)
      {
        CHECKRANGEIND(ind);

        pair<UINT, ULONG> loc;

        for(UINT i=1; i<=m_numSubPop; ++i)
        {
          if( m_subPopIndex[i] > ind )
          {
            loc.first = i-1;
            loc.second = ind - m_subPopIndex[i-1];
            break;
          }
        }
        return loc;
      }

      /// beginning index for subPopulation subPop
      /**
      \param subPop subpopulation index
      \result beginning index of this subpopulation
      \sa absIndIndex

      absIndIndex(index, subPop) = subPopBegin(subPop) + index
      */
      ULONG subPopBegin(UINT subPop) const
      {
        CHECKRANGESUBPOP(subPop);

        return m_subPopIndex[subPop];
      }

      /// ending index for subPopulation subPop
      /**
      \param subPop subpopulation index
      \result ending index of this subpopulation (not in this subpop)
      \sa absIndIndex
       \note as with all ...End functions, the returning index is out of the range
       so the actual range is [xxxBegin, xxxEnd). This agrees with all STL
       conventions.
      */
      ULONG subPopEnd(UINT subPop) const
      {
        CHECKRANGESUBPOP(subPop);

        return m_subPopIndex[subPop+1];
      }

      //@}
      /** @name itertors and accessers
        ways to access information, mainly various iterators.
      */
      //@{

      /// reference to individual ind
      /**
      \param ind absolute index of an individual
      \return reference to an individual
      */
      IndType& individual(ULONG ind)
      {
        CHECKRANGEIND(ind);
        DBG_ASSERT( ind < popSize(), IndexError,
          "Individual index " + toStr(ind) + " is out of range of 0 ~ " + toStr(popSize()-1));

        return m_inds[ind];
      }

      /// refernce to individual ind in subPopulation subPop
      /** \param ind individual index within subPop
          \param subPop subPopulation index
          \return reference to an individual
      */
      IndType& individual(ULONG ind, UINT subPop)
      {
        CHECKRANGESUBPOPMEMBER(ind, subPop);

        return m_inds[ subPopBegin(subPop) + ind];
      }

      /// CPPONLY individual iterator: without subPop info
      IndIterator indBegin()
      {
        return m_inds.begin();
      }

      /// CPPONLY individual iterator: without subPop info
      IndIterator indEnd()
      {
        return m_inds.end();
      }

      /// CPPONLY individual iterator: with subPop info.
      IndIterator indBegin(UINT subPop)
      {
        CHECKRANGESUBPOP(subPop);

        return m_inds.begin() + absIndIndex(0, subPop);
      }

      /// CPPONLY individual iterator: with subPop info.
      IndIterator indEnd(UINT subPop)
      {
        CHECKRANGESUBPOP(subPop);

        return m_inds.begin() + m_subPopIndex[subPop+1];
      }

      /// allele iterator that access a locus across all copies of chromosomes and individual
      /** CPPONLY
      \param locus
        allele access, given locus, return
       the first allele. ptr++ go the next one.
       default return the beginning of the first subPopulation,
       also the first of the whole population

      \note The order of alleles DOES NOT HAVE TO match the order of
      individuals. Only the boundary of subPopulations will be respected.
      Therefore, it is possible to access all alleles within an subPopulation
      through such iterators.
      */
      AlleleIterator alleleBegin(UINT locus)
      {
        CHECKRANGEABSLOCUS(locus);

        return AlleleIterator( &m_genotype[locus], totNumLoci());
      }

      /// CPPONLY allele iterator
      AlleleIterator alleleEnd(UINT locus)
      {
        CHECKRANGEABSLOCUS( locus);

        return AlleleIterator( &m_genotype[locus] + m_popGenoSize , totNumLoci());
      }

      ///  CPPONLY allele begin, for given subPopu
      AlleleIterator alleleBegin(UINT locus, UINT subPop)
      {
        CHECKRANGEABSLOCUS(locus);
        CHECKRANGESUBPOP(subPop);

        if(IndType::shallowCopiedFlagOn())
          adjustGenoPosition();

        return AlleleIterator( &m_genotype[ m_subPopIndex[subPop]*genoSize() +
          locus ], totNumLoci());
      }

      ///  CPPONLY allele iterator
      AlleleIterator alleleEnd( UINT locus, UINT subPop)
      {
        CHECKRANGEABSLOCUS(locus);
        CHECKRANGESUBPOP(subPop);

        if(IndType::shallowCopiedFlagOn())
          adjustGenoPosition();
        return AlleleIterator( &m_genotype[  m_subPopIndex[subPop+1]*genoSize() +
          locus ], totNumLoci());
      }

      ///  CPPONLY allele iterator, go through all allels one by one, without subPop info
      SimpleIterator begin()
      {
        return m_genotype.begin();
      }

      ///  CPPONLY allele iterator
      SimpleIterator end()
      {
        return m_genotype.end();
      }

      ///  CPPONLY allele iterator, go through all allels one by one in a subPopulation
      SimpleIterator begin(UINT subPop)
      {
        CHECKRANGESUBPOP(subPop);

        if(IndType::shallowCopiedFlagOn())
          adjustGenoPosition();

        return m_genotype.begin() + m_subPopIndex[subPop]*genoSize();
      }

      ///  CPPONLY allele iterator in a subPopulation.
      SimpleIterator end(UINT subPop)
      {
        CHECKRANGESUBPOP(subPop);

        if(IndType::shallowCopiedFlagOn())
          adjustGenoPosition();
        return m_genotype.begin() + m_subPopIndex[subPop+1]*genoSize();
      }

      ///  CPPONLY genoIterator --- beginning of individual ind.
      GenoIterator genoBegin(ULONG ind) const
      {
        CHECKRANGEIND(ind);
        return m_inds[ind].genoBegin();
      }

      ///  CPPONLY genoIterator -- end of individual ind.
      GenoIterator genoEnd(ULONG ind) const
      {
        CHECKRANGEIND(ind);
        return m_inds[ind].genoEnd();
      }

      ///  CPPONLY genoIterator --- beginning of individual ind.
      GenoIterator genoBegin(ULONG ind, UINT subPop) const
      {
        CHECKRANGESUBPOP(subPop);
        CHECKRANGESUBPOPMEMBER(ind, subPop);

        return m_inds[ subPopBegin(subPop) + ind].genoBegin();
      }

      ///  CPPONLY genoIterator -- end of individual ind.
      GenoIterator genoEnd(ULONG ind, UINT subPop) const
      {
        CHECKRANGESUBPOP(subPop);
        CHECKRANGESUBPOPMEMBER(ind, subPop);

        return m_inds[ subPopBegin(subPop) + ind].genoEnd();
      }

      /// get the whole genotype.
      /// individuals will be in order before exposing
      /// their genotypes.
      PyObject* arrGenotype()
      {
        if(IndType::shallowCopiedFlagOn())
          // adjust position. deep=true
          adjustGenoPosition(true);

        // false: directly expose values. Do not copy data over.
        return Allele_Vec_As_NumArray(m_popGenoSize, &m_genotype[0], false);
      }

      /// get the whole genotype.
      /// individuals will be in order before exposing
      /// their genotypes.
      PyObject* arrGenotype(UINT subPop)
      {
        CHECKRANGESUBPOP(subPop);
        if(IndType::shallowCopiedFlagOn())
          // adjust position. deep=true
          adjustGenoPosition(true);

        return Allele_Vec_As_NumArray(subPopSize(subPop)*genoSize(), &*begin(subPop), false);
      }

      //@}

      /** @name utility functions
       set subPopulation, save and load etc.
      */
      //@{

      /** brief return individual info in pop namespace
       */
      PyObject* exposeInfo()
      {
        // regardness of info type, treat it as int.
        vectori val(m_popSize);
        for(ULONG it=0; it < m_popSize; ++it)
          val[it] = static_cast<int>(individual(it).info());
        // this has to be changed according to info type.
        PyObject* var = setIntNumArrayVar("info", m_popSize, &val[0]);
        Py_INCREF(var);
        return var;
      }

      /** brief return individual affected status in pop namespace
       */
      PyObject* exposeAffectedness()
      {
        // regardness of info type, treat it as int.
        vectori val(m_popSize);
        for(ULONG it=0; it < m_popSize; ++it)
          val[it] = static_cast<int>(individual(it).affected());
        // this has to be changed according to info type.
        PyObject* var = setIntNumArrayVar("affected", m_popSize, &val[0]);
        Py_INCREF(var);
        return var;
      }

      // set info field of all individuals
      /**
      This function set info field of all individuals. Info can be used to
      sort individuals and set subpopulations. Therefore, the following code
      will do a migration:

      setIndInfo([ an_array_of_info_value ])

      setSubPopByIndInfo()

      \param info an array of info values, should have length of pop size
      \sa Individual::setInfo, Individual::info, info
      */
      void setIndInfo( std::vector<InfoType> info)
      {
        DBG_ASSERT( info.size() == m_popSize, ValueError,
          "Info should have the same length as pop size");

        for(ULONG it=0; it < m_popSize; ++it)
          individual(it).setInfo( info[it] );
      }

      /// set individual info with their subPopulation id.
      /**
      set individual info by subpop id.
      */
      void setIndInfoWithSubPopID()
      {
        for( UINT i=0, iEnd = numSubPop(); i < iEnd;  ++i)
          for( IndIterator it = indBegin(i), itEnd = indEnd(i); it != itEnd;  ++it)
            it -> setInfo(i);
      }

      /// adjust subPopulation according to individual info values
      /** assume individual has subpop index in their info value and put
      them into corresponding subpopulations.
      \param info optional info that will be used to set sub pop
      \note individual with negative info will be removed!
      \sa setInfo,
      */
      void setSubPopByIndInfo(vectori info=vectori())
      {
        if( !info.empty())
        {
          DBG_ASSERT( info.size() == m_popSize, ValueError,
            "Info should have the same length as pop size");
          for(ULONG it=0; it < m_popSize; ++it)
            individual(it).setInfo( info[it] );
        }

        DBG_DO(DBG_POPULATION, cout << "Sorting individuals."<< endl);
        // sort individuals first
        std::sort(indBegin(), indEnd());

        // sort individuals first
        // remove individuals with negative index.
        if( indBegin()->info() < 0 )
        {
          // popsize etc will be changed.
          ULONG newPopSize = m_popSize;
          IndIterator it=indBegin();
          for(; it != indEnd();  ++it)
          {
            if( it->info() < 0 )
              newPopSize -- ;
            else
              break;
          }
          // 'it' now point to the one with positive info()

          DBG_DO(DBG_POPULATION, cout << "New pop size" << newPopSize << endl);

          // allocate new genotype and inds
          ULONG newPopGenoSize = genoSize() * newPopSize;
          vectora newGenotype(newPopGenoSize);
          vector<IndType> newInds(newPopSize);

          DBG_ASSERT( indEnd()== newPopSize+it, SystemError,
            "Pointer misplaced. ");

          // assign genotype location and set structure information for individuals
          Allele * ptr = &newGenotype[0];
          UINT step = genoSize();
          for(ULONG i=0; i< newPopSize; ++i, ptr+=step, ++it)
          {
            newInds[i].setGenoStruIdx(genoStruIdx());
            newInds[i].setGenoPtr( ptr );
            newInds[i].copyFrom(*it);             // copy everything, with info value
          }

          // now, switch!
          m_genotype.swap(newGenotype);
          m_inds.swap(newInds);

          m_popSize = newPopSize;
          m_popGenoSize = newPopGenoSize;
        }

        if( m_inds.empty())
        {
          m_numSubPop = 1;
          m_subPopSize.resize(1,0);
          m_subPopIndex.resize(2);
        }
        else
        {
          // reset indices etc.
          m_numSubPop = static_cast<UINT>(m_inds.back().info())+1;
          m_subPopSize.resize(m_numSubPop);
          m_subPopIndex.resize(m_numSubPop+1);

          // check subpop size
          fill(m_subPopSize.begin(), m_subPopSize.end(), 0);
          for(IndIterator it = indBegin(), itEnd = indEnd(); it < itEnd;  ++it)
            m_subPopSize[ static_cast<UINT>(it->info()) ] ++;
        }
        /// rebuild index
        size_t i = 1;
        for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
          m_subPopIndex[i] = m_subPopIndex[i-1] + m_subPopSize[i - 1];
      }

      /// split population
      /** split subpopulation 'which' into subpopulations with size specified in subPops,
      optionally given subPopID.
      The subPOP ID of Other subpop will be kept. For example, if
      subpop 1 of 0 1 2 3 is split into three parts, the new subpop id will be
      0 (1 4 5) 2 3.
      \note subpop with negative id will be removed. So, you can shrink
      one subpop by split and set one of the new subpop with negative id.

      */
      void splitSubPop(UINT which, vectorlu sizes, vectori subPopID=vectori())
      {
        DBG_ASSERT( accumulate(sizes.begin(), sizes.end(), 0UL) == subPopSize(which),
          ValueError,
          "Sum of subpopulation sizes does not equal to the size of subpopulation to be splitted.");

        DBG_FAILIF( !subPopID.empty() && subPopID.size() != sizes.size(), ValueError,
          "If subPopID is given, it should have the same length as subPOP");

        if( sizes.size() == 1)
          return;

        // set initial info
        setIndInfoWithSubPopID();

        int spID;
        if(subPopID.empty())                      // starting sp number
          spID = which;
        else
          spID = subPopID[0];
        ULONG sz=0;                               // idx within subpop
        size_t newSPIdx=0;
        for(IndIterator ind = indBegin(which); ind != indEnd(which); ++ind)
        {
          if( sz == sizes[newSPIdx])
          {
            sz = 0;
            newSPIdx++;
            if(subPopID.empty())
              spID = numSubPop()+newSPIdx-1;
            else
              spID = subPopID[newSPIdx];
          }
          ind->setInfo(spID);
          sz++;
        }
        setSubPopByIndInfo();
      }

      /// split population
      /** split subpopulation 'which' into subpopulations with specified 'proportions',
      optionally given subPopID.
      \note subpop with negative id will be removed. So, you can shrink
      one subpop by split and set one of the new subpop with negative id.
      */
      void splitSubPopByProportion(UINT which, vectorf proportions, vectori subPopID=vectori())
      {
        DBG_ASSERT( fcmp_eq(accumulate(proportions.begin(), proportions.end(), 0.), 1), ValueError,
          "Proportions do not add up to one.");

        if( proportions.size() == 1)
          return;

        ULONG spSize = subPopSize(which);
        vectorlu subPop(proportions.size());
        for(size_t i=0; i< proportions.size()-1; ++i)
          subPop[i] = static_cast<ULONG>(spSize*proportions[i]);
        // to avoid round off problem, calculate the last subpopulation
        subPop[ subPop.size()-1] = spSize - accumulate( subPop.begin(), subPop.end()-1, 0L);
        splitSubPop(which, subPop, subPopID);
      }

      /** remove empty subpops, this will adjust subPOP ID of other subpops */
      void removeEmptySubPops()
      {
        // if remove empty subpops
        UINT newSPNum = m_numSubPop;
        vectorlu newSPSize;
        for(size_t sp=0; sp < m_numSubPop; ++sp)
        {
          if( m_subPopSize[sp] == 0 )
            newSPNum --;
          else
            newSPSize.push_back( m_subPopSize[sp]);
        }
        m_numSubPop = newSPNum;
        m_subPopSize.swap(newSPSize);
        m_subPopIndex.resize(m_numSubPop+1);
        /// rebuild index
        size_t i = 1;
        for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
          m_subPopIndex[i] = m_subPopIndex[i-1] + m_subPopSize[i - 1];
      }

      /**  remove subpop, adjust subpop numbers so that there will be no 'empty'
      subpops left */
      void removeSubPops(const vectoru& subPops=vectoru(), bool removeEmptySubPops=false)
      {
        setIndInfoWithSubPopID();
        int shift=0;
        for( size_t sp = 0; sp < m_numSubPop; ++sp)
        {
          if( find( subPops.begin(), subPops.end(), sp) != subPops.end())
          {
            shift++;
            for(IndIterator ind = indBegin(sp); ind != indEnd(sp); ++ind)
              ind->setInfo(-1);                   // remove
          }
          // other subpop shift left
          else
          {
            for(IndIterator ind = indBegin(sp); ind != indEnd(sp); ++ind)
              ind->setInfo(sp-shift);             // shift left
          }
        }
        setSubPopByIndInfo();
        if(removeEmptySubPops)
          this->removeEmptySubPops();
      }

      /**  remove subpop, adjust subpop numbers so that there will be no 'empty'
      subpops left */
      void removeIndividuals(const vectoru& inds=vectoru())
      {
        setIndInfoWithSubPopID();
        for(size_t i = 0; i < inds.size(); ++i)
          individual(i).setInfo(-1);              // remove
        setSubPopByIndInfo();
      }

      /// merge population
      /** merge subpopulations, subpop id will be the ID of the first in array subPops */
      void mergeSubPops(vectoru subPops=vectoru(), bool removeEmptySubPops=false)
      {
        // set initial info
        setIndInfoWithSubPopID();

        // merge all subpopulations
        if(subPops.empty())
        {
          // [ popSize() ]
          vectorlu sz(1, popSize());
          setSubPopStru(sz, false);
          return;
        }

        UINT id = subPops[0];
        for(UINT sp=0; sp < numSubPop(); ++sp)
        {
          if( find( subPops.begin(), subPops.end(), sp) != subPops.end())
            for(IndIterator ind = indBegin(sp); ind != indEnd(sp); ++ind)
              ind->setInfo(id);
        }
        setSubPopByIndInfo();

        if( removeEmptySubPops)
          this->removeEmptySubPops();
      }

      /// \brief reorder subpopulations
      /**
      \param order new order of the subpopulations. For examples, 3 2 0 1
      means subpop3, subpop2, subpop0, subpop1 will be the new layout.
      \param rank you can also specify new rank for each subpop. For example, 3,2,0,1
      means the original subpopulations will have new ID 3,2,0,1. To achive order 3,2,0,1.
      the rank should be 1 0 2 3.
      */
      void reorderSubPops(const vectoru& order=vectoru(), const vectoru& rank=vectoru(),
        bool removeEmptySubPops=false)
      {

        DBG_FAILIF( order.empty() && rank.empty(), ValueError,
          "Please specify one of order or rank.");

        DBG_FAILIF( !order.empty() && !rank.empty(), ValueError,
          "You can specify only one of order or rank.");

        if(removeEmptySubPops)
          this->removeEmptySubPops();

        if( ( !order.empty() && order.size() != m_numSubPop )
          || ( ! rank.empty() && rank.size() != m_numSubPop))
          cout << "Warning: Given order or rank does not have the length of number of subpop." << endl;

        if( !order.empty())
        {
          // alow order[i] > numSubPop(). In a special case, I have last empty subpop...
          for(size_t i=0; i < order.size(); ++i)
          {
            if(order[i] >= numSubPop())
              continue;
            for(IndIterator ind = indBegin(order[i]); ind != indEnd(order[i]); ++ind)
              ind->setInfo(i);
          }
        }
        else
        {
          for(size_t i=0; i < rank.size(); ++i)
          {
            if(i >= numSubPop())
              continue;
            for(IndIterator ind = indBegin(i); ind != indEnd(i); ++ind)
              ind->setInfo(rank[i]);
          }
        }
        // reset ...
        setSubPopByIndInfo();
      }

      /** form a new population according to info */
      Population<Ind>& newPopByIndInfo(bool keepAncestralPops=true, vectori info=vectori(), bool removeEmptySubPops=false)
      {
        // copy the population over (info is also copied)
        Population<Ind>& pop = this->clone(keepAncestralPops);
        // and shrink them
        pop.setSubPopByIndInfo(info);
        if( removeEmptySubPops)
          pop.removeEmptySubPops();
        return pop;
      }

      void removeLoci( const vectoru& remove=vectoru(), const vectoru& keep=vectoru())
      {
        DBG_FAILIF( !keep.empty() && ! remove.empty(), ValueError,
          "Please specify one and only one of keep or remove.");

        if( keep.empty() && remove.empty() )
          return;

        vectoru loci;
        if( !keep.empty())
          loci = keep;
        else
        {
          for(size_t loc = 0; loc < this->totNumLoci(); ++loc)
            // if not removed
            if( find(remove.begin(), remove.end(), loc) == remove.end())
              loci.push_back(loc);
        }

#ifndef OPTIMIZED
        for(size_t i=0; i<loci.size(); ++i)
        {
          DBG_FAILIF( loci[i] >= this->totNumLoci(), ValueError,
            "Given loci " + toStr(loci[i]) + " exceed max number of loci." );
          DBG_FAILIF( i > 0 && loci[i] <= loci[i-1], ValueError,
            "Given loci should be in order.");
        }
#endif
        // adjust order before doing anything
        UINT newTotNumLoci = loci.size();
        UINT oldTotNumLoci = this->totNumLoci();

        // first, new genotype
        // get a GenoStructure with parameters. GenoStructure may be shared by some populations
        // a whole set of functions ploidy() etc in GenoStruTriat can be used after this step.
        vectoru newNumLoci;
        vectorf newLociDist;
        vectorstr newLociNames;
        UINT lastCh = -1;
        for(vectoru::iterator loc = loci.begin();
          loc != loci.end(); ++loc)
        {
          UINT ch = this->chromLocusPair(*loc).first;
          if( newNumLoci.empty() || lastCh != ch )
          {
            newNumLoci.push_back(1);
            lastCh = ch;
          }
          else
            newNumLoci.back()++;
          newLociDist.push_back( this->locusPos(*loc));
          newLociNames.push_back( this->locusName(*loc));
        }

        // new geno structure is in effective now!
        this->setGenoStructure(this->ploidy(), newNumLoci, this->sexChrom(), newLociDist,
          this->alleleNames(), newLociNames, this->maxAllele() );
        // prepare data
        //
        // keep m_popSize;
        // keep m_numSubPop;
        // keep m_subPopSize;
        // keep m_subPopIndex;

        // genotype
        // allocate new genotype and inds
        ULONG newPopGenoSize = genoSize() * m_popSize;
        vectora newGenotype(newPopGenoSize);
        // keep newInds();

        // copy data over
        Allele * ptr = &newGenotype[0];
        UINT pEnd = this->ploidy();
        for(ULONG i=0; i< m_popSize; ++i)
        {
          // set new geno structure
          m_inds[i].setGenoStruIdx(genoStruIdx());
          Allele * oldPtr = m_inds[i].genoPtr();
          // new genotype
          m_inds[i].setGenoPtr( ptr );
          // copy each chromosome
          for(UINT p=0; p < pEnd; ++p)
          {
            for(vectoru::iterator loc = loci.begin();
              loc != loci.end(); ++loc)
            {
              *(ptr++) = oldPtr[*loc];
            }
            oldPtr += oldTotNumLoci;              // next ploidy
          }
        }
        m_popGenoSize = newPopGenoSize;
        m_genotype.swap(newGenotype);

        // ancestral populations?
        for(size_t ap=0; ap < m_ancestralPops.size(); ++ap)
        {
          popData& p = m_ancestralPops[ap];
          // set pointers
          vector<IndType>& inds = p.m_inds;
          ULONG ps = inds.size();
          vectora newGenotype(ps*pEnd*newTotNumLoci);
          ptr = &(newGenotype[0]);

          for(ULONG i=0; i< ps; ++i)
          {
            // set new geno structure
            inds[i].setGenoStruIdx(genoStruIdx());
            Allele * oldPtr = inds[i].genoPtr();
            // new genotype
            inds[i].setGenoPtr( ptr );
            // copy each chromosome
            for(UINT p=0; p < pEnd; ++p)
            {
              for(vectoru::iterator loc = loci.begin();
                loc != loci.end(); ++loc)
              {
                *(ptr++) = oldPtr[*loc];
              }
              oldPtr += oldTotNumLoci;            // next ploidy
            }
          }
          p.m_genotype.swap(newGenotype);
        }                                         // all ancestral
      }

      /** get a new population with selected loci */
      Population<Ind>& newPopWithPartialLoci(
        const vectoru& remove=vectoru(),
        const vectoru& keep=vectoru())
      {
        // copy the population over (info is also copied)
        Population<Ind>* pop = new Population<Ind>(*this);
        pop->removeLoci(remove, keep);
        return *pop;
      }

      // swap in rhs. (usually a scratch population
      // if fixRhs is false, leave rhs in broken status
      // if force is true, push current population
      // regardless of m_ancestray settings.
      void pushAndDiscard(Population<Ind>& rhs, bool force=false)
      {
        // time consuming!
        DBG_ASSERT( rhs.genoStru() == genoStru(), ValueError,
          "Passed population has different genotypic structure");

        DBG_ASSERT( &(m_genotype[0]) != &(rhs.m_genotype[0]), ValueError,
          "Passed population is a reference of current population, swapPop failed." );

        // front -1 pop, -2 pop, .... end
        //
        if( !force && m_ancestralDepth > 0
          && ancestralDepth() == static_cast<size_t>(m_ancestralDepth) )
          m_ancestralPops.pop_back();

        // save current population
        if( force || m_ancestralDepth != 0 )
        {
          // add a empty popData
          m_ancestralPops.push_front(popData());
          // get its reference
          popData& pd = m_ancestralPops.front();
          // swap with real data
          // current population may *not* be in order
          pd.m_subPopSize.swap(m_subPopSize);
          // store starting geno ptr,
          // if m_genotype is re-allocated, reset pointers
          // in m_inds
#ifndef OPTIMIZED
          pd.m_startingGenoPtr = &m_genotype[0];
#endif
          pd.m_genotype.swap(m_genotype);
          pd.m_inds.swap(m_inds);
        }

        // then swap out data
#ifndef OPTIMIZED
        Allele* rhsStartingGenoPtr = &(rhs.m_genotype[0]);
        Allele* lhsStartingGenoPtr = &(m_genotype[0]);
#endif
        m_popSize = rhs.m_popSize;
        m_numSubPop = rhs.m_numSubPop;
        m_subPopSize.swap( rhs.m_subPopSize);
        m_popGenoSize = rhs.m_popGenoSize;
        m_subPopIndex.swap( rhs.m_subPopIndex);
        m_genotype.swap( rhs.m_genotype);
        m_inds.swap(rhs.m_inds);
#ifndef OPTIMIZED
        DBG_FAILIF( rhsStartingGenoPtr != &(m_genotype[0]),
          SystemError, "Starting genoptr has been changed.");
        DBG_FAILIF( lhsStartingGenoPtr != &(rhs.m_genotype[0]),
          SystemError, "Starting genoptr has been changed.");
#endif
        // current population should be working well
        // (with all datamember copied form rhs
        // rhs may not be working well since m_genotype etc
        // may be from ancestral pops
        if( rhs.m_popSize != rhs.m_inds.size())
        {
          // keep size if pop size is OK.
          // remove all supopulation structure of rhs
          rhs.m_popSize = rhs.m_inds.size();
          rhs.m_numSubPop = 1;
          rhs.m_subPopSize.resize(1, rhs.m_popSize);
          rhs.m_popGenoSize = genoSize()*rhs.m_popSize;
          rhs.m_subPopIndex.resize(2,0);
          rhs.m_subPopIndex[1] = rhs.m_popSize;
          // no need to set genoPtr or genoStru()
        }
      }

      UINT ancestralDepth()
      {
        return m_ancestralPops.size();
      }

      /// set ancestral depth, can be -1
      void setAncestralDepth(int depth)
      {
        // just to make sure.
        useAncestralPop(0);
        //
        if(depth >=0 && m_ancestralPops.size() > static_cast<size_t>(depth))
        {
          int numRemove = m_ancestralPops.size() - depth;
          while( numRemove-- > 0 )
            m_ancestralPops.pop_back();
        }
        DBG_ASSERT( depth <0 || m_ancestralPops.size() <= static_cast<size_t>(depth), SystemError,
          "Failed to change ancestral Depth");

        m_ancestralDepth = depth;
      }

      int ancestralPop()
      {
        return m_curAncestralPop;
      }

      // idx = 0 (current, do nothing) -1, -2, ..
      //
      void useAncestralPop(int idx)
      {
        if( idx == m_curAncestralPop)
          return;

        DBG_DO(DBG_POPULATION, cout << "Use ancestralPop: " << idx <<
          "Curidx: " <<  m_curAncestralPop << endl);

        if( idx == 0 || m_curAncestralPop != 0)   // recover pop.
        {
          popData& pd = m_ancestralPops[ m_curAncestralPop-1 ];
          pd.m_subPopSize.swap(m_subPopSize);
          pd.m_genotype.swap(m_genotype);
          pd.m_inds.swap(m_inds);
          m_curAncestralPop = 0;
#ifndef OPTIMIZED
          DBG_FAILIF( pd.m_startingGenoPtr != &m_genotype[0],
            SystemError, "Starting genoptr has been changed.");
          pd.m_startingGenoPtr = &(pd.m_genotype[0]);
#endif
          if( idx == 0)
          {                                       // restore key paraemeters from data
            m_popSize = m_inds.size();
            m_popGenoSize = genoSize() * m_popSize;
            m_numSubPop = m_subPopSize.size();
            m_subPopIndex.resize( m_numSubPop + 1);
            // build subPop index
            UINT i = 1;
            for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
              m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];
            return;
          }
        }

        // now m_curAncestralPop is zero.
        DBG_ASSERT( idx >0 && static_cast<size_t>(idx-1) < m_ancestralPops.size(),
          ValueError, "Ancestry population " + toStr(idx) + " does not exist.");

        // now idx should be at least 1
        m_curAncestralPop = idx;
        // swap  1 ==> 0, 2 ==> 1

        popData& pd = m_ancestralPops[ m_curAncestralPop -1];
        pd.m_subPopSize.swap(m_subPopSize);
        pd.m_genotype.swap(m_genotype);
        pd.m_inds.swap(m_inds);
#ifndef OPTIMIZED
        pd.m_startingGenoPtr = &(pd.m_genotype[0]);
#endif
        // use pd
        m_popSize = m_inds.size();
        m_popGenoSize = genoSize() * m_popSize;
        m_numSubPop = m_subPopSize.size();
        m_subPopIndex.resize( m_numSubPop + 1);
        UINT i = 1;
        for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
          m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];

      }

      /// compare two populations
      bool equalTo(const Population<Ind>& rhs)
      {
        return(
          genoStru() == rhs.genoStru() &&
          m_subPopSize == rhs.m_subPopSize &&
          m_inds == rhs.m_inds );
      }

      //@}

      /// CPPONLY
      /// some iterators requires that genotype information is within
      /// each subPopulation. We need to adjust genotypic info to
      /// obey this.
      void adjustGenoPosition(bool deep=false);

      /// save population to a file
      /**
      \param filename save to filename
      \param format format to save. Can be one of 'text', 'bin', 'xml'

      The default format is 'text' but the output is not suppored to be read.
      'bin' has smaller size and should be used for large populations.
      'xml' format is most readable and should be used when you would
      like to convert simuPOP populations to other formats.
      \sa global function loadPopulation
      */
      void savePopulation(const string& filename, const string& format="auto")
      {
        ofstream ofs(filename.c_str());

        if(!ofs)
          throw ValueError("Can not open file " + filename );

        if( format == "text" || (format == "auto" && filename.substr(filename.size()-4,4) == ".txt" ))
        {
          boost::archive::text_oarchive oa(ofs);
          oa << *this;
        }
#ifndef __NO_XML_SUPPORT__
        else if (format == "xml" || (format == "auto" && filename.substr(filename.size()-4,4) == ".xml" ) )
        {
          boost::archive::xml_oarchive oa(ofs);
          oa << boost::serialization::make_nvp("population",*this);
        }
#endif
        else if (format == "bin" ||  (format == "auto" && filename.substr(filename.size()-4,4) == ".bin" ))
        {
          boost::archive::binary_oarchive oa(ofs);
          oa << *this;
        }
        else
        {
          ofs.close();
          throw ValueError("Wrong format type. Use one of text, xml, bin or appropriate extension txt, xml or bin");
        }
        ofs.close();

      }

      /// CPPONLY load population from a file
      /**
      \param filename load from filename
      \param format format to load. Can be one of "text", "bin", "xml". It should match
        the format used to save the population.

      \sa savePopulation
      */
      void loadPopulation(const string& filename, const string& format="auto")
      {

        ifstream ifs(filename.c_str());

        if(!ifs)
          throw ValueError("Can not open file " + filename );

        if( format == "text" || (format == "auto" && filename.substr(filename.size()-4,4) == ".txt" ))
        {
          boost::archive::text_iarchive ia(ifs);
          ia >> *this;
        }
#ifndef __NO_XML_SUPPORT__
        else if (format == "xml" ||  (format == "auto" && filename.substr(filename.size()-4,4) == ".xml" ))
        {
          boost::archive::xml_iarchive ia(ifs);
          ia >> boost::serialization::make_nvp("population",*this);
        }
#endif
        else if (format == "bin" || (format == "auto" && filename.substr(filename.size()-4,4) == ".bin" ) )
        {
          boost::archive::binary_iarchive ia(ifs);
          ia >> *this;
        }
        else
        {
          ifs.close();
          throw ValueError("Wrong format type. Use one of text, xml, bin or appropriate extension txt, xml or bin");
        }
        ifs.close();
      }

    public:

      int rep()
      {
        return m_rep;
      }

      /// CPPONLY  set rep number
      void setRep(int rep, bool setVar=true)
      {
        m_rep = rep;
        if(setVar)
          m_vars.setIntVar("rep", rep);
      }

      int grp()
      {
        return m_grp;
      }

      /// CPPONLY
      void setGrp(int grp, bool setVar=true)
      {
        m_grp = grp;
        if(setVar)
          m_vars.setIntVar("grp", grp);
      }

      ULONG gen()
      {
        return m_gen;
      }

      /// CPPONLY
      void setGen(ULONG gen, bool setVar=true)
      {
        m_gen = gen;
        if(setVar)
          m_vars.setIntVar("gen", gen);
      }

      /// return variables of this population
      /// if subPop is given, return dictionary for
      /// specified subpopulation.
      PyObject* vars(int subPop=-1)
      {
        if(subPop < 0)
        {
          Py_INCREF(m_vars.dict());
          return m_vars.dict();
        }
        else
        {
          DBG_ASSERT( static_cast<UINT>(subPop) < numSubPop() ,
            IndexError, "Subpop index out of range of 0 ~ " + toStr(numSubPop()-1) );

          DBG_ASSERT( hasVar("subPop"), ValueError,
            "subPop statistics does not exist yet.");

          PyObject* spObj = getVar("subPop");
          spObj = PyList_GetItem(spObj, subPop);

          DBG_ASSERT(spObj != NULL, SystemError,
            "Something is wrong about the length of subPop list. ");

          Py_INCREF( spObj);
          return spObj;
        }
      }

      /// CPPONLY
      /// The same as vars(), but without increasing
      /// reference count.
      PyObject* dict(int subPop=-1)
      {
        if(subPop < 0)
          return m_vars.dict();
        else
        {
          DBG_ASSERT( static_cast<UINT>(subPop) < numSubPop() ,
            IndexError, "Subpop index out of range of 0 ~ " + toStr(numSubPop()-1) );

          DBG_ASSERT( hasVar("subPop"), ValueError,
            "subPop statistics does not exist yet.");

          PyObject* spObj =  getVar("subPop");
          spObj = PyList_GetItem(spObj, subPop);

          DBG_ASSERT(spObj != NULL, SystemError,
            "Something is wrong about the length of subPop list. ");

          return spObj;
        }
      }

      /// CPPONLY
      void setDict(PyObject* dict)
      {
        DBG_ASSERT(dict != NULL, SystemError, "Dictionary is empty");
        m_vars.setDict(dict);
      }

      ///
      bool hasVar(const string& name)
      {
        return m_vars.hasVar(name);
      }

      /// CPPNLY
      void removeVar(const string& name)
      {
        m_vars.removeVar(name);
      }

#ifndef OPTIMIZED
      void checkRefCount()
      {
        m_vars.checkRefCount();
      }
#endif

      /// CPPONLY
      PyObject* setBoolVar(const string& name, const bool val)
      {
        return  m_vars.setBoolVar(name, val);
      }

      /// CPPONLY
      PyObject* setIntVar(const string& name, const int val)
      {
        return m_vars.setIntVar(name, val);
      }

      /// CPPONLY
      PyObject* setDoubleVar(const string& name, const double val)
      {
        return m_vars.setDoubleVar(name, val);
      }

      /// CPPONLY
      PyObject* setStringVar(const string& name, const string& val)
      {
        return m_vars.setStringVar(name, val);
      }

      /// CPPONLY
      PyObject* setStrDictVar(const string& name, const strDict& val)
      {
        return m_vars.setStrDictVar(name, val);
      }

      /// CPPONLY
      PyObject* setIntDictVar(const string& name, const intDict& val)
      {
        return m_vars.setIntDictVar(name, val);
      }

      /// CPPONLY
      PyObject* setVar(const string& name, PyObject* val)
      {
        return m_vars.setVar(name, val);
      }

      /// CPPONLY
      PyObject* setDoubleNumArrayVar(const string& name, int dim, double*buf, bool copyOver=true)
      {
        return m_vars.setDoubleNumArrayVar(name, dim, buf, copyOver);
      }

      /// CPPONLY
      PyObject* setIntNumArrayVar(const string& name, int dim, int*buf, bool copyOver=true)
      {
        return m_vars.setIntNumArrayVar(name, dim, buf, copyOver);
      }

      /// CPPONLY
      PyObject* getVar(const string&name, bool nameError=true)
      {
        return m_vars.getVar(name, nameError);
      }

      /// CPPONLY
      bool getVarAsBool(const string& name, bool nameError=true)
      {
        return m_vars.getVarAsBool(name, nameError);
      }

      /// CPPONLY
      int getVarAsInt(const string& name, bool nameError=true)
      {
        return m_vars.getVarAsInt(name, nameError);
      }

      /// CPPONLY
      double getVarAsDouble(const string& name, bool nameError=true)
      {
        return m_vars.getVarAsDouble(name, nameError);
      }

      /// CPPONLY
      string getVarAsString(const string& name, bool nameError=true)
      {
        return m_vars.getVarAsString(name, nameError);
      }

      /// CPPONLY
      strDict getVarAsStrDict(const string& name, bool nameError=true)
      {
        return m_vars.getVarAsStrDict(name, nameError);
      }

      /// CPPONLY
      intDict getVarAsIntDict(const string& name, bool nameError=true)
      {
        return m_vars.getVarAsIntDict(name, nameError);
      }

      /// CPPONLY
      int getVarAsDoubleNumArray(const string& name, double* & buf,  bool nameError=true)
      {

        return m_vars.getVarAsDoubleNumArray(name, buf, nameError);
      }

      /// CPPONLY
      int getVarAsIntNumArray(const string& name, int* & buf,  bool nameError=true)
      {
        return m_vars.getVarAsIntNumArray(name, buf, nameError);
      }

      /// CPPONLY
      string varsAsString() const
      {
        return m_vars.asString();
      }

      /// CPPONLY
      void varsFromString(const string& vars)
      {
        return m_vars.fromString(vars);
      }

      /// evaluate python statment/expressions
      /**
        this function evaluate python expressions
        and return as string representing the result
          */
      string evaluate(const string& expr="", const string& stmts="")
      {
        return Expression(expr, stmts, m_vars.dict() ).valueAsString();
      }

      ///
      void execute(const string& stmts="")
      {
        Expression("", stmts, m_vars.dict() ).evaluate();
      }

    private:

      friend class boost::serialization::access;

      template<class Archive>
        void save(Archive &ar, const UINT version) const
      {
        // deep adjustment: everyone in order
        const_cast<Population*>(this)->adjustGenoPosition(true);

        ar & make_nvp("libraryMaxAllele", MaxAllele);
        ar & make_nvp("geno_structure", this->genoStru());
        ar & make_nvp("subPop_sizes", m_subPopSize);
        ar & make_nvp("genotype", m_genotype);
        ar & make_nvp("individuals", m_inds);
        ar & make_nvp("ancestry", m_ancestralDepth);
        size_t sz = m_ancestralPops.size();
        ar & make_nvp("numOfAncestralPops", sz);
        for(size_t i=0; i< m_ancestralPops.size(); ++i)
        {
          const_cast<Population*>(this)->useAncestralPop(i+1);
          // need to make sure ancestral pop also in order
          const_cast<Population*>(this)->adjustGenoPosition(true);
          ar & make_nvp("subPop_sizes", m_subPopSize);
          ar & make_nvp("genotype", m_genotype);
          ar & make_nvp("individuals", m_inds);
        }
        const_cast<Population*>(this)->useAncestralPop(0);

        // save shared variables as string.
        // note that many format are not supported.
        string vars = varsAsString();
        ar & make_nvp("vars", vars);
      }

      template<class Archive>
        void load(Archive &ar, const UINT version)
      {
        ULONG ma;
        ar & make_nvp("libraryMaxAllele", ma);

        if( ma > MaxAllele)
          cout << "Warning: The population is saved in long allele library. \n"
            << "Unless all alleles are less than "  << MaxAllele
            << ", you should use the long allele library. (c.f. simuOpt.setOptions()\n";

        GenoStructure  stru;
        ar & make_nvp("geno_structure", stru);
        ar & make_nvp("subPop_sizes", m_subPopSize);
        ar & make_nvp("genotype", m_genotype);
        ar & make_nvp("individuals", m_inds);

        // set genostructure, check duplication
        // we can not use setGenoStruIdx since stru may be new.
        this->setGenoStructure(stru);

        m_numSubPop = m_subPopSize.size();
        m_popSize = accumulate(m_subPopSize.begin(), m_subPopSize.end(), 0L);

        if( m_popSize != m_inds.size() )
          throw ValueError("Number of individuals does not match population size.");

        m_popGenoSize = totNumLoci() * GenoStruTrait::ploidy()
          * m_popSize;

        m_subPopIndex.resize(m_numSubPop + 1);
        UINT i = 1;
        for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
          m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];

        // assign genotype location and set structure information for individuals
        Allele * ptr = &(m_genotype[0]);
        UINT step = genoSize();
        for(ULONG i=0; i< m_popSize; ++i, ptr+=step)
        {
          m_inds[i].setGenoStruIdx(genoStruIdx());
          m_inds[i].setGenoPtr( ptr );
        }
        m_ancestralDepth = 0;
        m_ancestralPops.clear();

        // ancestry populations
        ar & make_nvp("ancestry", m_ancestralDepth);
        size_t na;
        ar & make_nvp("numOfAncestralPops", na);
        for(size_t ap=0; ap< na; ++ap)
        {
          popData pd;
          ar & make_nvp("subPop_sizes", pd.m_subPopSize);
          ar & make_nvp("genotype", pd.m_genotype);
          ar & make_nvp("individuals", pd.m_inds);
          // set pointer after copy this thing again (push_back)
          m_ancestralPops.push_back(pd);
          // now set pointers
          popData& p = m_ancestralPops.back();
          // set pointers
          vector<IndType>& inds = p.m_inds;
          ULONG ps = inds.size();
          ptr = &(p.m_genotype[0]);

          for(ULONG i=0; i< ps; ++i, ptr+=step)
          {
            inds[i].setGenoPtr( ptr );
            // set new genoStructure
            inds[i].setGenoStruIdx(genoStruIdx());
            // fresh copy so clear shallowcopied flag.
            inds[i].setShallowCopied(false);
          }
        }

        // load vars from string
        string vars;
        ar & make_nvp("vars", vars);
        varsFromString(vars);
      }

      BOOST_SERIALIZATION_SPLIT_MEMBER();

    private:
      /// population size: number of individual
      ULONG m_popSize;

      /// number of subPopulations
      UINT m_numSubPop;

      /// size of each subPopulation
      vectorlu m_subPopSize;

      ///  size of genotypic information of the whole poopulation
      LONG m_popGenoSize;

      /// index to subPop \todo change to vectorl
      vectorlu m_subPopIndex;

      /// pool of genotypic information
      vectora m_genotype;

      /// individuals.
      vector<IndType> m_inds;

      int m_ancestralDepth;

      /// shared variables for this population
      SharedVariables m_vars;

      /// store previous populations
      /// need to store: subPopSize, genotype and m_inds
      struct popData
      {
#ifndef OPTIMIZED
        Allele* m_startingGenoPtr;
#endif
        vectorlu m_subPopSize;
        vectora m_genotype;
        vector<IndType> m_inds;
      };

      std::deque<popData> m_ancestralPops;

      /// curent replicate, group number
      int m_rep, m_grp;

      /// generation
      ULONG m_gen;

      /// current ancestral depth
      int m_curAncestralPop;

  };

  /// CPPONLY
  template< class Ind>
    void Population<Ind>::adjustGenoPosition(bool deep)
  {

    // everyone in strict order
    if(deep)
    {
      vectora tmpGenotype(m_popGenoSize);
      vectora::iterator it = tmpGenotype.begin();

      for(IndIterator ind=indBegin(), indEd=indEnd(); ind!=indEd; ++ind)
      {
        copy(ind->genoBegin(), ind->genoEnd(), it);
        ind->setGenoPtr(&*it);
        it+=this->totNumLoci()*this->ploidy();
        ind->setShallowCopied(false);
      }
      // discard original genotype
      tmpGenotype.swap(m_genotype);
      // set geno pointer
      IndType::clearShallowCopiedFlag();
      return;
    }

    /// find out how many individuals are shallow copied.
    vectorl scIndex(0);
    ULONG j, k=0, jEnd;

    for(UINT sp=0, spEd = numSubPop(); sp < spEd;  sp++)
    {
      Allele * spBegin = &*m_genotype.begin() + m_subPopIndex[sp]*genoSize();
      Allele * spEnd   = &*m_genotype.begin() + m_subPopIndex[sp+1]*genoSize();
      for( j=0, jEnd = subPopSize(sp); j < jEnd;  j++)
      {
        if( m_inds[k].shallowCopied() )
        {
          if( &*genoBegin(k) < spBegin
            || &*genoEnd(k) > spEnd )
            /// record individual index and genoPtr
            scIndex.push_back(k);
          else
            m_inds[k].setShallowCopied(false);
        }
        k++;
      }
    }

    if( scIndex.empty()) return;

    /// to further save time, deal with a special case that there are
    /// only two shallowCopied individuals
    if( scIndex.size() == 2 )
    {
      // swap!
      Allele* tmp = m_inds[ scIndex[0] ].genoPtr();
      m_inds[ scIndex[0] ].setGenoPtr( m_inds[ scIndex[1] ].genoPtr() );
      m_inds[ scIndex[1] ].setGenoPtr( tmp);

      Allele tmp1;
      for(UINT a =0; a < genoSize(); a++)
      {
        tmp1 =  m_inds[ scIndex[0] ].allele(a);
        m_inds[ scIndex[0] ].setAllele(a, m_inds[ scIndex[1] ].allele(a) );
        m_inds[ scIndex[1] ].setAllele(a, tmp1);
      }

      m_inds[ scIndex[0] ].setShallowCopied(false);
      m_inds[ scIndex[1] ].setShallowCopied(false);
      IndType::clearShallowCopiedFlag();
      return;
    }

    /// save genotypic info
    vectora scGeno(scIndex.size() * totNumLoci() * ploidy());
    vector<Allele*> scPtr( scIndex.size() );

    size_t i, iEnd;

    for(i=0, iEnd = scIndex.size(); i < iEnd;  i++)
    {
      scPtr[i] = m_inds[ scIndex[i]].genoPtr();
      copy( genoBegin(scIndex[i]), genoEnd(scIndex[i]), scGeno.begin() + i* genoSize());
    }

    DBG_DO(DBG_POPULATION, cout << "Shallow copied" << scIndex << endl);

    /// sort the pointers!
    sort( scPtr.begin(), scPtr.end());

    /// copy back.
    for(i=0, iEnd =scIndex.size(); i < iEnd;  i++)
    {
      m_inds[ scIndex[i] ].setGenoPtr( scPtr[i]);
      copy( scGeno.begin() + i*  genoSize(), scGeno.begin() + (i+1)*  genoSize(),
        genoBegin( scIndex[i] ));
      m_inds[ scIndex[i] ].setShallowCopied(false);
    }

    IndType::clearShallowCopiedFlag();
  }

}
#endif
