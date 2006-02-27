/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$
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

#ifndef _INITIALIZER_H
#define _INITIALIZER_H
/**
\file
\brief head file of class initializer:public Operator
*/
#include "utility.h"
#include "operator.h"

#include <numeric>
using std::string;
using std::accumulate;

#include "Python.h"

namespace simuPOP
{

  /**
  \brief  initialize alleles at the start of generation.

  @author Bo Peng
  */

  class initializer: public Operator
  {
    public:
      /// constructor. default to be always active.
      initializer( const vectoru& subPop=vectoru(),
        intMatrix indRange=intMatrix(),
        const vectoru& atLoci = vectoru(),
        int atPloidy = -1,
        double maleFreq=0.5, const vectori& sex = vectori(),
        int stage=PreMating, int begin=0, int end=-1, int step=1,
        vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL)
        :Operator("","", stage, begin, end, step, at, rep, grp),
        m_subPop(subPop), m_indRange(indRange),
        m_atLoci(atLoci), m_atPloidy(atPloidy),
        m_maleFreq(maleFreq), m_sex(sex)
      {
        for(size_t i = 0; i < m_indRange.size(); ++i)
        {
          // allow for singleton
          if( m_indRange[i].size() == 1)
            m_indRange[i].push_back( m_indRange[i][0]+1);

          if( m_indRange[i].size() != 2 || m_indRange[i][0] > m_indRange[i][1] )
            throw ValueError("Expecting a range.");
        }
      }

      /// destructor
      virtual ~initializer()
      {
      }

      /// this function is very important
      virtual Operator* clone() const
      {
        return new initializer(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::initializer>";
      }

      void setRanges(population& pop)
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
              m_ranges[i].resize(2, pop.subPopBegin( m_subPop[i]));
              m_ranges[i][1] = pop.subPopEnd(m_subPop[i]);
            }
          }
        }
      }

      void initSexIter()
      {
        if(! m_sex.empty())
          m_sexItr = m_sex.begin();
      }

      Sex nextSex()
      {
        DBG_ASSERT(*m_sexItr == int(Male) || *m_sexItr == int(Female),
          ValueError, "sex must be array of Male or Female. "
          + toStr(*m_sexItr) + " given.");
        Sex s = *m_sexItr++ == 1?Male:Female;
        if(m_sexItr == m_sex.end())
          m_sexItr = m_sex.begin();
        return s;
      }

    protected:

      /// applicable subpop
      vectoru m_subPop;

      /// ranges
      intMatrix m_indRange;

      /// init loci
      vectoru m_atLoci;

      /// at which ploidy, -1 means all
      int m_atPloidy;

      /// sex frequency
      double m_maleFreq;

      /// specify sex
      vectori m_sex;

      /// iterator to sex
      vectori::iterator m_sexItr;

      /// populaiton specific range
      intMatrix m_ranges;

  };

  /// initialize genotype by allele frequency and sex by male frequency

  class initByFreq: public initializer
  {
    public:
      /** \brief randomly assign alleles according to allele frequency

      This operator randomly assign alleles according to given allele frequency.
      Allele frequencies can differ by subpop. Sex is also assigned randomly.

      \param alleleFreq an array of allele frequencies. Must add up to 1; or
         a matrix of allele frequencies, each row corresponse
         a subpopulation.
      \param subPop an array of applicable subpopulations. default to all
      \param indRange a [begin, end] pair of range of individuals; or
        an array of [begin, end] pairs.
      \param identicalInds whether or not make individual genotype identical
      in all subpopulation. If true, this operator will randomly generate genotype for
      an individual and spread it to the whole subpopulation.
      \param atLoci a vector of loci indices. If empty, apply to all loci
      \param atPloidy initialize which copy of chromosomes. Default to all.
      \param maleFreq male frequency. Default to 0.5.
      \param sex an arry of sex [Male, Female, Male]... for individuals.
      The length of sex will not be checked. If length of sex is shorter than
      number of individuals, sex will be reused from the beginning.
      \param stages is set to PreMating. Other parameters please see help(baseOperator.__init__)
      */
      initByFreq( const matrix& alleleFreq=matrix(),
        bool identicalInds=false,  const vectoru& subPop=vectoru(),
        intMatrix indRange = intMatrix(),
        const vectoru& atLoci=vectoru(), int atPloidy=-1,
        double maleFreq=0.5, const vectori& sex = vectori(),
        int stage=PreMating, int begin=0, int end=1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        : initializer(subPop, indRange, atLoci,
        atPloidy, maleFreq, sex,
        stage, begin, end, step, at, rep, grp),
        m_alleleFreq(alleleFreq), m_identicalInds(identicalInds)
      {

        DBG_FAILIF(  m_alleleFreq.empty(),
          IndexError, "Should specify one of alleleFreq, alleleFreqs");

        for(size_t i=0; i< m_alleleFreq.size(); ++i)
          if( fcmp_ne(accumulate(m_alleleFreq[i].begin(), m_alleleFreq[i].end(), 0.), 1.0))
            throw ValueError("Allele frequencies should add up to one.");
      }

      ~initByFreq()
      {
      }

      /// this function is very important
      virtual Operator* clone() const
      {
        return new initByFreq(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::initByFreq>";
      }

      bool apply(population& pop)
      {
        /// initialize m_ranges
        setRanges(pop);

        this->initSexIter();

        DBG_FAILIF( m_alleleFreq.size() > 1 && m_alleleFreq.size() != this->m_ranges.size(),
          ValueError, "Ranges and values should have the same length");

        for(size_t rg = 0; rg < this->m_ranges.size(); ++rg)
        {
          vectorf& alleleFreq = m_alleleFreq.size() == 1 ? m_alleleFreq[0] : m_alleleFreq[rg];

          ULONG left = this->m_ranges[rg][0], right = this->m_ranges[rg][1];

          DBG_FAILIF( left > pop.popSize() || right > pop.popSize() || left > right ,
            ValueError, "Invaid range boundary: " + toStr(left) + " - " + toStr(right-1));

          // Weightedsampler ws(rng(), incFreq);
          Weightedsampler ws(rng(), alleleFreq);

          DBG_ASSERT( fcmp_eq(std::accumulate(alleleFreq.begin(), alleleFreq.end(), 0.), 1),
            SystemError, "Allele frequecies shoudl add up to one.");

          pop.adjustGenoPosition(true);
          if( m_identicalInds )
          {
            if( this->m_atLoci.empty())           // to all loci
            {
              if( this->m_atPloidy==-1)           // all chromosomes
              {
                ws.get(pop.indGenoBegin(left), pop.indGenoEnd(left), StartingAllele);

                for(ULONG ind=left+1; ind != right; ++ind)
                  copy(pop.indGenoBegin(left), pop.indGenoEnd(left), pop.indGenoBegin(ind));
              }
              else                                // only initialize one set of chromosome
              {
                ws.get(pop.ind(left).genoBegin(this->m_atPloidy),
                  pop.ind(left).genoEnd(this->m_atPloidy), StartingAllele);

                for(ULONG ind=left+1; ind != right; ++ind)
                  copy(pop.ind(left).genoBegin(this->m_atPloidy),
                    pop.ind(left).genoEnd(this->m_atPloidy),
                    pop.ind(ind).genoBegin(this->m_atPloidy));
              }
            }
            else
            {
              // initislize locus by locus
              for(vectoru::iterator locus=this->m_atLoci.begin(); locus != this->m_atLoci.end(); ++locus)
              {
                if( this->m_atPloidy == -1 )      // all chromosomes
                {
                  ws.get( pop.alleleBegin(*locus) + left * pop.ploidy(),
                    pop.alleleBegin(*locus) + (left+1)*pop.ploidy(), StartingAllele );

                  for(ULONG ind=left+1; ind != right; ++ind)
                    copy(pop.alleleBegin(*locus) + left*pop.ploidy(),
                      pop.alleleBegin(*locus) + (left+1)*pop.ploidy(),
                      pop.alleleBegin(*locus) + ind*pop.ploidy());
                }
                else                              // only one of the copies (do it one by one?)
                {
                  UINT a = ws.get() + StartingAllele;
                  for(ULONG ind=left; ind != right; ++ind)
                    pop.ind(ind).setAllele(a, *locus, this->m_atPloidy);
                }
              }
            }
          }
          else                                    // not idential individuals
          {
            if( this->m_atLoci.empty())           // at all loci
            {
              if( this->m_atPloidy == -1)
                ws.get( pop.genoBegin()+left*pop.genoSize(),
                  pop.genoBegin()+right*pop.genoSize(), StartingAllele);
              else                                // for only one ploidy
              {
                for(ULONG ind=left; ind != right; ++ind)
                  ws.get( pop.ind(ind).genoBegin(this->m_atPloidy),
                    pop.ind(ind).genoEnd(this->m_atPloidy), StartingAllele);
              }
            }
            else                                  // at certain loci
            {
              if( this->m_atPloidy == -1)
              {
                for(vectoru::iterator locus=this->m_atLoci.begin(); locus != this->m_atLoci.end(); ++locus)
                  ws.get( pop.alleleBegin(*locus) + left * pop.ploidy(),
                    pop.alleleBegin(*locus) + right * pop.ploidy(), StartingAllele );
              }
              else                                // for only one ploidy
              {
                for(ULONG ind=left; ind != right; ++ind)
                  for(vectoru::iterator locus=this->m_atLoci.begin(); locus != this->m_atLoci.end(); ++locus)
                    pop.ind(ind).setAllele(ws.get()+StartingAllele, *locus, this->m_atPloidy);
              }
            }
          }
          // initialize sex
          // call randUnif once for each individual
          // (initialize allele need to call randUnif for each locus
          if(this->m_sex.empty())
          {
            for ( ULONG ind = left; ind != right; ++ind)
            {
              if( rng().randUniform01() < this->m_maleFreq )
                pop.ind(ind).setSex( Male );
              else
                pop.ind(ind).setSex( Female );
            }
          }
          else                                    // use provided sex array
          {
            for ( ULONG ind = left; ind != right; ++ind)
            {
              pop.ind(ind).setSex( this->nextSex());
            }
          }                                       // set sex
        }                                         // range
        return true;
      }

    private:

      /// allele frequencies (assume all loci are the same for a subPop
      matrix m_alleleFreq;

      ///
      bool m_identicalInds;

  };

  /// initialize genotype by value and then copy to all individuals

  class initByValue:public initializer
  {
    public:
      /*** \brief initialize populations by given alleles. Every individual will have the same genotype.

      This operator assign given alleles to specified individuals. The parameter combinations should be

      value - subPop / indRange / indRanges : individual in subPops or in range/ranges will
      be assigned genotype 'value'
         subPop /indRanges: subPop or indRanges should have the same length as values. Each item
      of values will be assigned to each subPop or indRange.

      \param value an array of genotypes of one individual, having the same length as
      the length of atLoci() or atLoci()*ploidy() or pop.genoSize() (whole genotype) or totNumLoci()
      (one copy of chromosome); or an array of array of genotypes of one individual. Should have length one or
      equal to subpop or ranges or proportion.
      \param atLoci a vector of loci indices. If empty, apply to all loci
      \param atPloidy initialize which copy of chromosomes. Default to all.
      \param subPop an array of applicable subpopulations. If values are given,
      should have equal length to values.
      \param indRange a [begin, end] pair of range of individuals; or
      an array of [begin, end] pairs.
      \param proportions an array of percentages for each item in values.
      \param maleFreq male frequency
      \param sex an arry of sex [Male, Female, Male]... for individuals.
      The length of sex will not be checked. If length of sex is shorter than
      number of individuals, sex will be reused from the beginning.
      \param stages is set to PreMating. Other parameters please see help(baseOperator.__init__)
      */
      initByValue( intMatrix value=intMatrix(),
        vectoru atLoci=vectoru(), int atPloidy=-1,
        vectoru subPop=vectoru(), intMatrix indRange=intMatrix(),
        const vectorf& proportions = vectorf(),
        double maleFreq=0.5, const vectori& sex = vectori(),
        int stage=PreMating, int begin=0, int end=1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        : initializer(subPop, indRange, atLoci, atPloidy, maleFreq, sex,
        stage, begin, end, step, at, rep, grp),
        m_value(value), m_proportion(proportions)
      {
        DBG_FAILIF( maleFreq < 0 || maleFreq > 1 ,
          IndexError, "male frequency in the population should be in the range of [0,1]");

        DBG_FAILIF( m_value.empty() , ValueError,
          "Please specify an array of alleles in the order of chrom_1...chrom_n for all copies of chromosomes");

        DBG_FAILIF( !m_proportion.empty() && m_proportion.size() != m_value.size(), ValueError,
          "If proportions are given, its length should match that of values.");

        DBG_FAILIF( !m_proportion.empty() && fcmp_ne(accumulate(m_proportion.begin(), m_proportion.end(), 0.0), 1),
          ValueError, "Proportion should add up to one.");
      }

      ~initByValue(){}

      /// this function is very important
      virtual Operator* clone() const
      {
        return new initByValue(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::initByValue>";
      }

      bool apply(population& pop)
      {

        this->initSexIter();

        for(size_t i = 0; i < m_value.size(); ++i)
        {
          // extend m_value
          if(m_proportion.empty())
          {
            if( this->m_atLoci.empty() && m_value[i].size() == pop.totNumLoci())
            {
              if( this->m_atPloidy==-1)
              {
                m_value[i].resize( pop.totNumLoci() * pop.ploidy());
                for(UINT p=1; p < pop.ploidy(); ++p)
                  copy( m_value[i].begin(), m_value[i].begin() + pop.totNumLoci(),
                    m_value[i].begin() + pop.totNumLoci()*p);
              }
            }
            if( ! this->m_atLoci.empty() && m_value[i].size() == this->m_atLoci.size() )
            {
              if( this->m_atPloidy==-1)
              {
                m_value[i].resize( this->m_atLoci.size() * pop.ploidy());
                for(UINT p=1; p < pop.ploidy(); ++p)
                  copy( m_value[i].begin(), m_value[i].begin() + this->m_atLoci.size(),
                    m_value[i].begin() + this->m_atLoci.size()*p);
              }
            }                                     // case 2
          }                                       // no proportion
        }                                         // check m_value[i]

        DBG_DO(DBG_INITIALIZER, cout << "Size of src is " << m_value[0].size() << endl);

#ifndef OPTIMIZED
        UINT gSz= m_value[0].size();
        for(size_t v=1; v< m_value.size(); ++v)
          if( m_value[v].size() != gSz)
            throw ValueError("Given values should have the same length (either one copy of chromosomes or the whole genotype.");
#endif

        setRanges(pop);

        DBG_FAILIF( m_proportion.empty() && m_value.size() > 1
          && m_value.size() != this->m_ranges.size(),
          ValueError, "Ranges and values should have the same length");

        if( m_proportion.empty() )
        {
          // for each range
          for(size_t rg = 0; rg < this->m_ranges.size(); ++rg)
          {
            // we have left and right
            ULONG left = this->m_ranges[rg][0], right = this->m_ranges[rg][1];

            DBG_FAILIF( left > pop.popSize() || right > pop.popSize() || left > right ,
              ValueError, "Invaid m_ranges boundary: " + toStr(left) + " - " + toStr(right));

            vectori& src = m_value.size() > 1 ? m_value[rg] : m_value[0];
            size_t srcSz = src.size(), lociSz=this->m_atLoci.size(),
              totNumLoci = pop.totNumLoci();

            for (ULONG ind = left; ind < right; ++ind)
            {
              if(this->m_atLoci.empty())
              {
                if(this->m_atPloidy==-1)          // all copies of chromosome
                {
                  DBG_ASSERT(src.size() == pop.genoSize(), ValueError,
                    "Length of value does not match geno size");
                  copy(src.begin(), src.end(), pop.indGenoBegin(ind));
                }
                else                              // one of the copied.
                {                                 /// fixme: check length of src?
                  DBG_ASSERT( src.size() == pop.totNumLoci(), ValueError,
                    "Ploidy is specified but the length of alleles do not match length of chromosome. val size: " );
                  copy(src.begin(), src.end(), pop.ind(ind).genoBegin(this->m_atPloidy));
                }
              }
              else                                // with m_loci
              {
                if(this->m_atPloidy==-1)          // all copied of chromosome
                {
                  DBG_ASSERT( src.size() == this->m_atLoci.size() ||
                    (src.size() == this->m_atLoci.size() * pop.ploidy() ),
                    ValueError, "Length of value does not atLoci size");
                  for(size_t loc = 0; loc != srcSz ;++loc)
                    *(pop.indGenoBegin(ind) + this->m_atLoci[loc%lociSz] + loc/lociSz*totNumLoci ) = src[loc];
                }
                else                              // one of the copies.
                {
                  DBG_ASSERT( src.size() == this->m_atLoci.size(), ValueError,
                    "Ploidy is specified but the length of alleles do not match length of given allele array.");
                  for(size_t loc = 0; loc != srcSz ;++loc)
                    *(pop.ind(ind).genoBegin(this->m_atPloidy) +
                    this->m_atLoci[loc%lociSz] + loc/lociSz*totNumLoci ) = src[loc];
                }
              }
              if( this->m_sex.empty())
              {
                if( rng().randUniform01() < this->m_maleFreq )
                  pop.ind(ind).setSex( Male );
                else
                  pop.ind(ind).setSex( Female );
              }
              else
              {
                pop.ind(ind).setSex( this->nextSex() );
              }                                   // set sex
            }
          }
        }
        else                                      // use proportion.
        {
          for(size_t rg = 0; rg < this->m_ranges.size(); ++rg)
          {
            ULONG left = this->m_ranges[rg][0], right = this->m_ranges[rg][1];

            DBG_FAILIF( left > pop.popSize() || right > pop.popSize() || left > right ,
              ValueError, "Invaid m_ranges boundary: " + toStr(left) + " - " + toStr(right));

            size_t srcSz = m_value[0].size(), lociSz=this->m_atLoci.size(), totNumLoci = pop.totNumLoci();
            Weightedsampler ws(rng(), m_proportion);

            for (ULONG ind = left; ind < right; ++ind)
            {
              if(this->m_atLoci.empty())          // whole chromosome or geno
              {
                if( srcSz == totNumLoci)          // by ploidy
                {
                  if(this->m_atPloidy==-1)
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
                      pop.indGenoBegin(ind)+this->m_atPloidy*totNumLoci);
                  }
                }
                else                              // whole geno
                {
                  UINT idx = ws.get();
                  if(this->m_atPloidy==-1)
                    copy(m_value[idx].begin(), m_value[idx].end(), pop.indGenoBegin(ind));
                  else                            // only one copy of chromosome
                    copy(m_value[idx].begin(), m_value[idx].end(),
                      pop.ind(ind).genoBegin(this->m_atPloidy));
                }
              }
              else                                /// atLoci is in effect
              {
                if( srcSz == lociSz )             // one by one
                {
                  if(this->m_atPloidy==-1)
                  {
                    for(UINT p = 0; p<pop.ploidy(); ++p)
                    {
                      UINT idx = ws.get();
                      for(size_t loc = 0; loc != srcSz ;++loc)
                        *(pop.indGenoBegin(ind) + p*totNumLoci + this->m_atLoci[loc]) = m_value[idx][loc];
                    }
                  }
                  else
                  {
                    UINT idx = ws.get();
                    for(size_t loc = 0; loc != srcSz ;++loc)
                      *(pop.indGenoBegin(ind) + this->m_atPloidy*totNumLoci +
                      this->m_atLoci[loc]) = m_value[idx][loc];
                  }

                }
                else                              // who geno (at loci .. though)
                {
                  if(this->m_atPloidy==-1)
                  {
                    UINT idx = ws.get();
                    for(size_t loc = 0; loc != srcSz ;++loc)
                      *(pop.indGenoBegin(ind) + this->m_atLoci[loc%lociSz] +
                      loc/lociSz*totNumLoci ) = m_value[idx][loc];
                  }
                  else
                  {
                    UINT idx = ws.get();
                    for(size_t loc = 0; loc != srcSz ;++loc)
                      *(pop.ind(ind).genoBegin(this->m_atPloidy)
                      + this->m_atLoci[loc%lociSz] +
                      loc/lociSz*totNumLoci ) = m_value[idx][loc];
                  }
                }
              }
              if( this->m_sex.empty())
              {
                if( rng().randUniform01() < this->m_maleFreq )
                  pop.ind(ind).setSex( Male );
                else
                  pop.ind(ind).setSex( Female );
              }
              else
              {
                pop.ind(ind).setSex( this->nextSex() );
              }                                   // set sex
            }
          }

        }

        return true;
      }

    private:
      /// allele frequencies (assume all loci are the same
      intMatrix m_value;

      /// if assign randomly
      vectorf m_proportion;
  };

  /// initialize genotype by value and then copy to all individuals

  class spread:public Operator
  {
    public:
      // copy genotype of ind to all individuals in subPop.
      spread(ULONG ind, vectoru subPop=vectoru(),
        int stage=PreMating, int begin=0, int end=1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        : Operator("","", stage, begin, end, step, at, rep, grp),
        m_ind(ind), m_subPop(subPop)
      {
      }

      ~spread(){}

      /// this function is very important
      virtual Operator* clone() const
      {
        return new spread(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::spread genotype>";
      }

      bool apply(population& pop)
      {
        std::pair<UINT, ULONG> p = pop.subPopIndPair(m_ind);

        if(m_subPop.empty())
          m_subPop.resize(1, p.first);

        GenoIterator srcBegin = pop.indGenoBegin(m_ind),
          srcEnd=pop.indGenoEnd(m_ind);

        for(vectoru::iterator sp=m_subPop.begin(); sp != m_subPop.end(); ++sp)
        {
          for(ULONG i = pop.subPopBegin(*sp); i < pop.subPopEnd(*sp); ++i)
            if( i != m_ind)
              copy(srcBegin, srcEnd, pop.indGenoBegin(i));
        }

        return true;
      }

    private:
      ULONG m_ind;
      vectoru m_subPop;

  };

  class pyInit:public initializer
  {
    public:
      /** \brief initialize populations using given user function.

      User of this operator must supply a Python function with parameter (index, ploidy, subpop).
      This operator will loop through all individual in each subpop and call this function
      to initialize populations.

      The arrange of parameters allows different initialization scheme for each subpop.

      \param func a python function with parameter (index, ploidy, subpop) index is the allele
      index (0 ~ totNumLoci()-1), ploidy (index to copy of chromosomes), subpop (subpop number).
      The return value of this function should be a integer.
      \param atLoci a vector of loci indices. If empty, apply to all loci
      \param atPloidy initialize which copy of chromosomes. Default to all.
      \param stage is et to PreMating. Other parameters please refer to help(baseOperator.__init__)
      */
      pyInit(PyObject * func,  vectoru subPop=vectoru(),
        vectoru atLoci=vectoru(), int atPloidy=-1,
        intMatrix indRange=intMatrix(),
        double maleFreq=0.5, const vectori& sex = vectori(),
        int stage=PreMating, int begin=0, int end=1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL)
        : initializer(subPop, indRange, atLoci, atPloidy, maleFreq, sex,
        stage, begin, end, step, at, rep, grp)
      {
        DBG_FAILIF( maleFreq < 0 || maleFreq > 1 ,
          IndexError, "male frequency in the population should be in the range of [0,1]");
        DBG_ASSERT( PyCallable_Check(func),
          ValueError, "Func is not a Python function");

        Py_XINCREF(func);
        m_func = func;
      }

      ~pyInit()
      {
        if( m_func != NULL)
          Py_DECREF(m_func);
      }

      /// CPPONLY
      pyInit(const pyInit& rhs):initializer(rhs), m_func(rhs.m_func)
      {
        if( m_func != NULL)
          Py_INCREF(m_func);
      }

      /// this function is very important
      virtual Operator* clone() const
      {
        return new pyInit(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::pyInit>";
      }

      bool apply(population& pop)
      {
        this->initSexIter();

        for(UINT sp=0, numSP=pop.numSubPop(); sp < numSP; ++sp)
        {
          for(ULONG it=0, itEnd=pop.subPopSize(sp); it<itEnd; ++it)
          {
            for(UINT p=0, pEnd=pop.ploidy(); p<pEnd; ++p)
            {
              for(UINT al = 0, alEnd=pop.totNumLoci(); al < alEnd; ++al)
              {
                int resInt;
                PyCallFunc3( m_func, "(iii)", al, p, sp, resInt, PyObj_As_Int);
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

    private:
      /// the python function with parameter (ind, ploidy, subpop)
      PyObject* m_func;

  };

}
#endif
