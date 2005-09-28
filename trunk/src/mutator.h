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

#ifndef _MUTATOR_H
#define _MUTATOR_H
/**
\file
\brief head file of class Mutator:public Operator
*/
#include "operator.h"

/// for hybrid mutator
#include "Python.h"

namespace simuPOP
{
  /** \brief mutator class.

  Do not use this class directly. It just provide interface for real mutators.

  Every mutator can specify rate (equal rate) or rates (different rate for different
  loci) and a vector of applicable loci (default to all but should have the same
  length with rates if rates have length greater than one).

  max allele can be specified as well but more parameter, if needed, should
  be implemented by individual mutator classes.

  Number of possible allelic states: Most theoretical studies assume an infinite
  number of allelic states to avoid any homoplasy. If it facilitates analysis,
  this is however extremely unrealistic.

  @author Bo Peng
  */
  template<class Pop>
    class Mutator: public Operator<Pop>
  {
    public:
      /** \brief create a mutator
      All mutators have the following common parameters. However, the actual meaning
      of these parameters may vary according to different model. Check the manual
      for details!!! (help(kamMutator) for example.)

      \param rate single rate for all applicable loci (atLoci). Will be ignored
      if rates is specified; or it can be an array of rates, the same length as atLoci.
      \param atLoci a vector of loci index. Can be ignored only when single
      rate is specified. Default to all loci.
      \param maxAllele max allowable allele. Interpreted by each sub mutaor class. Default to pop.maxAllele().
      */
      Mutator( vectorf rate=vectorf(),
        vectori atLoci=vectori(), UINT maxAllele=0,
        string output=">", string outputExpr="", 
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t")
        :Operator<Pop>(output, outputExpr, stage, begin, end, step, at, rep, grp, sep), 
        m_rate(rate), m_maxAllele(maxAllele), m_atLoci(atLoci),
        m_bt(rng()), m_initialized(false)
      {
        if( m_rate.empty() )
          throw SystemError("You should specify at least rate ( (0,1] ) or rates (a sequence of rate.)");

        if( rate.size() > 1 && atLoci.empty())
          throw SystemError("If you use variable rates, you should specify atLoci for each of the rate.");

        if( rate.size() > 1 && !atLoci.empty() && rate.size() != atLoci.size() )
          throw SystemError("If both rates and atLoci are specified, they should have the same length.");

        DBG_ASSERT( maxAllele <= MaxAllele, ValueError,
          "The maximum allele number exceeds " + toStr(MaxAllele)
          + ". \nIf you need longer allele size, please use simuPOP_la libraries.");
      }

      /// destructor
      virtual ~Mutator()
      {
      }

      /// this function is very important
      virtual Operator<Pop>* clone() const
      {
        return new Mutator<Pop>(*this);
      }

      /// return mutation rate
      vectorf rate()
      {
        return m_rate;
      }

      /// set an equal mutation rate
      void setRate(const double rate, const vectori atLoci = vectori() )
      {
        m_rate.resize(1);
        m_rate[0] = rate;
        if( ! atLoci.empty())
          m_atLoci = atLoci;

        m_initialized = false;
      }

      /// set an array of rates
      void setRate(const vectorf rate, const vectori atLoci = vectori() )
      {
        if( rate.size() != 1 && rate.size() != atLoci.size() )
          throw ValueError("If you specify more than one rate values, you should also specify corresponding applicable loci");

        m_rate = rate;
        if( ! atLoci.empty())
          m_atLoci = atLoci;

        m_initialized = false;
      }

      /// return max allowable allele number
      UINT maxAllele()
      {
        return m_maxAllele;
      }

      ///
      void setMaxAllele(UINT maxAllele)
      {
        m_maxAllele = maxAllele;
      }

      /// how to mutate a single allele.
      /// this is usually the only function that need to be defined by the subclasses.
      virtual void mutate(Allele* allele)
      {
        throw SystemError("You are not supposed to call this base mutator funciton.");
      };

      /// apply!
      virtual bool apply(Pop& pop)
      {
        if( !m_initialized || m_bt.size() != pop.ploidy() * pop.popSize())
        {
          initialize(pop);
          DBG_DO(DBG_MUTATOR, cout << "Mutate at loci " << m_atLoci <<
            " at rate " << m_rate << endl);
        }

#ifndef OPTIMIZED
        vectori mutated(m_atLoci.size());
        fill(mutated.begin(), mutated.end(), 0);
#endif

        DBG_DO(DBG_MUTATOR, cout <<"Mutate replicate " << pop.rep() << endl);
        DBG_ASSERT(pop.rep() == pop.getVarAsInt("rep"), SystemError,
          "Replicate number mismatch");

        m_bt.doTrial();

        // mutate each mutable locus
        for( size_t i=0, iEnd=m_atLoci.size(); i < iEnd; ++i)
        {
          const BitSet& succ = m_bt.succ(i);

          int locus = m_atLoci[i];

          BitSet::size_type pos = succ.find_first();
          if( pos != BitSet::npos)
          {
            do
            {
#ifndef OPTIMIZED
              Allele* ptr = (pop.alleleBegin( locus ) + pos ).ptr();
              if( ptr < &*pop.begin() || ptr >= &*pop.end() )
                throw SystemError("Mutation happens on an allele out of range.") ;
              DBG_DO(DBG_MUTATOR, cout << "Mutate locus " << locus
                << " of individual " << (pos/pop.ploidy()) << " from " << int(*ptr) );
              mutate(ptr);
              DBG_DO(DBG_MUTATOR, cout << " to " << int(*ptr) << endl);
              mutated[i]++;
#else
              mutate( (pop.alleleBegin( locus ) + pos ).ptr());
#endif
            }while( (pos = succ.find_next(pos)) != BitSet::npos );
          }                                       // succ.any
        }                                         // each applicable loci

        DBG_DO( DBG_MUTATOR, cout << "Mutated alleles: " << mutated << endl);

        return true;
      }

    private:

      /// initialize bernulli trial according to pop size etc
      virtual void initialize(Pop& pop)
      {
        if( m_maxAllele == 0 )
          m_maxAllele = pop.maxAllele();
        else if ( m_maxAllele > 0 && m_maxAllele > pop.maxAllele() )
          throw ValueError("maxAllele exceeds population max allele.");

        DBG_DO(DBG_MUTATOR, cout << "initialize mutator" << endl);

        // deal with applicable loci
        if(m_atLoci.empty() )
        {
          // all loci
          m_atLoci.resize(pop.totNumLoci() );
          for(UINT i=0, iEnd = pop.totNumLoci(); i < iEnd;  ++i)
            m_atLoci[i] = i;
        }

        /// all use the same rate
        if( m_rate.size() < m_atLoci.size() )
        {
          m_rate.resize( m_atLoci.size());
          fill(m_rate.begin()+1, m_rate.end(), m_rate[0]);
        }

        m_bt.setParameter(m_rate, pop.ploidy() * pop.popSize());

#ifndef OPTIMIZED
        for(size_t i=0; i<m_rate.size(); ++i)
          if( fcmp_lt( m_rate[i], 0.) || fcmp_gt( m_rate[i], 1.) )
            throw ValueError("Migration rate should be between [0,1], given " + toStr(m_rate[i]));
#endif

        m_initialized = true;
      }

    private:
      /// mutation rates
      vectorf m_rate;

      /// maxAllele
      UINT m_maxAllele;

      /// applicable loci.
      vectori m_atLoci;

      /// bernulli trials. bitSet mutation results.
      BernulliTrials m_bt;

      /// initialized? the first apply() call will trigger an initialization process.
      bool m_initialized;

  };

  /// K-Allele Model mutator
  /**
  Under this model, there are K (here refers as maxAllele) possible allele states, and any
  allele has a constant probability (rate/(K-1)) of mutation towards any of the K-1 allelic states.

  \note the theoretical mutation rate is rates/(K-1)  towards any of the K-1 allelic states.
  So rates is actually the probability to mutate!

  \sa Crow & Kimura 1970
  */

  template<class Pop>
    class KAMMutator: public Mutator<Pop>
  {
    public:
      /// K-Allele Model mutator
      /**
      \param rate  mutation rate. It is 'probability to mutate'. The actual
         mutation rate to any of the other K-1 allelic states are rates/(K-1)!
      \param states allelic states. Default to is 1 - maxAllele . Given states should be in order!
      \param atLoci and other parameters: refer to help(mutator), help(baseOperator.__init__)
      */
      KAMMutator(vectorf rate=vectorf(), vectori atLoci=vectori(),
        UINT maxAllele=0, vectora states=vectora(),
        string output=">", string outputExpr="", 
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t")
        :Mutator<Pop>( rate, atLoci, maxAllele,
        output, outputExpr, stage, begin, end, step, at, rep, grp, sep),
        m_states(states)
      {
        if( !m_states.empty())
        {
          DBG_FAILIF( m_states.size() == 1,
            ValueError, "Mutation should have more than one states.");

          DBG_FAILIF( m_states[0] < 1, ValueError, "Allele states must be >= 1");

          for( size_t i=1; i< m_states.size(); ++i)
          {
            DBG_ASSERT( m_states[i] > m_states[i-1], ValueError,
              "Given allelic states should be in order.");
          }
        }
      }

      ~KAMMutator(){}

      /// mutate to a state other than current state with equal probability
      virtual void mutate(Allele* allele)
      {
        if(m_states.empty())
        {
          Allele new_allele = rng().randInt(this->maxAllele()-1)+1;
          if(new_allele >= *allele)
            *allele = new_allele+1;
          else
            *allele = new_allele;
        }
        else
        {
          DBG_FAILIF( find(m_states.begin(), m_states.end(), *allele) == m_states.end(),
            ValueError, "Allelic state is not in one of the specified states.");

          size_t idx = rng().randInt(m_states.size()-1);
          if( m_states[idx] >= *allele)
            *allele = m_states[idx+1];
          else
            *allele = m_states[idx];
        }
      }

      /// this function is very important
      virtual Operator<Pop>* clone() const
      {
        return new KAMMutator<Pop>(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::k-allele model mutator K=" +
          toStr(m_states.empty()?this->maxAllele():m_states.size()) + ">" ;
      }

    private:
      vectora m_states;

  };

  /// stepwise mutation model.
  /**
  Stepwise mutation model (SMM) assumes that alleles are represented by integer values
  and that a mutation either increases or decreases the allele value by one.
  For variable number tandem repeats loci (VNTR), the allele value is generally
  taken as the number of tandem repeats in the DNA sequence.

  \sa Kimura & Ohta 1978
  */
  template<class Pop>
    class SMMMutator: public Mutator<Pop>
  {
    public:
      ///
      /**
      The stepwise mutation model (SMM) is developed for allozymes. It  provides better description
       for these kinds of evolutionary processes.

      \param rate: mutation rate
      \param incProb probability to increase allele state. Default to 1
      \param atLoci and other parameters: refer to help(mutator), help(baseOperator.__init__)

      */
      SMMMutator(vectorf rate=vectorf(), vectori atLoci=vectori(),
        UINT maxAllele=0, double incProb=0.5,
        string output=">", string outputExpr="", 
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t")
        :Mutator<Pop>( rate, atLoci, maxAllele,
        output, outputExpr, stage, begin, end, step, at, rep, grp, sep),
        m_incProb(incProb)
      {
        DBG_ASSERT( fcmp_ge( incProb, 0.) && fcmp_le( incProb, 1.),
          ValueError, "Inc probability should be between [0,1], given " + toStr(incProb));
      }

      ~SMMMutator(){}

      virtual void mutate(Allele* allele)
      {
        if( rng().randUniform01() < m_incProb)
        {
          if( *allele < this->maxAllele() )
            ++*allele;
        }
        else
        {
          if( *allele > 1 )
            --*allele;
        }
      }

      /// this function is very important
      virtual Operator<Pop>* clone() const
      {
        return new SMMMutator<Pop>(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::step-wise mutation model mutator>" ;
      }

    private:
      /// probability to increase allele state
      double m_incProb;
  };

  /// stepwise mutation model.
  /**
  Generalized Stepwise mutation model (GSM) assumes that alleles are represented by integer values
  and that a mutation either increases or decreases the allele value by a random value.

  \sa Kimura & Ohta 1978
  */
  template<class Pop>
    class GSMMutator: public Mutator<Pop>
  {
    public:
      ///
      /**
      The generalized stepwise mutation model (GMM) is developed for allozymes.
      It  provides better description
       for these kinds of evolutionary processes.

      \param rate: mutation rate
      \param incProb probability to increase allele state. Default to 1
      \param atLoci and other parameters: refer to help(mutator), help(baseOperator.__init__)
      \param func return number of steps. no parameter
      */
      GSMMutator( vectorf rate=vectorf(), vectori atLoci=vectori(),
        UINT maxAllele=0, double incProb=0.5, double p=0, PyObject* func=NULL,
        string output=">", string outputExpr="", 
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t")
        :Mutator<Pop>( rate, atLoci, maxAllele,
        output, outputExpr, stage, begin, end, step, at, rep, grp, sep),
        m_incProb(incProb), m_p(p), m_func(func)
      {
        DBG_ASSERT( fcmp_ge( incProb, 0.) && fcmp_le( incProb, 1.),
          ValueError, "Inc probability should be between [0,1], given " + toStr(incProb));

        if( func != NULL)                         // use this function
        {
          DBG_ASSERT( PyCallable_Check(func),
            ValueError, "Func is not a Python function");

          Py_XINCREF(func);
          m_func = func;
        }
        else
        {
          DBG_ASSERT( fcmp_ge( p, 0.) && fcmp_le( p, 1.),
            ValueError, "Parameter p of a geometric distribution should be between [0,1], given " + toStr(m_p));

        }
      }

      ~GSMMutator(){}

      /// this function is very important
      virtual Operator<Pop>* clone() const
      {
        return new GSMMutator<Pop>(*this);
      }

      virtual void mutate(Allele* allele)
      {
        int step;

        if( m_func == NULL)                       // use a geometric distribution.
          step = rng().randGeometric(m_p);
        else
        {
          PyObject* arglist = Py_BuildValue("()");
          PyObject* result = PyEval_CallObject(m_func, arglist);
          Py_DECREF(arglist);
          if( result == NULL)
          {
            PyErr_Print();
            throw ValueError("Function call failed.");
          }

          PyObj_As_Int(result, step);
          Py_DECREF(result);
        }

        DBG_DO(DBG_MUTATOR, cout << "step is " << step << endl);

        if( rng().randUniform01() < m_incProb)
        {
          if( static_cast<UINT>(*allele + step) < this->maxAllele() )
            *allele += step;
          else
            *allele = this->maxAllele();
        }
        else
        {
          if( *allele - step > 1 )
            *allele -= step;
          else
            *allele = 1;
        }
      }

      virtual string __repr__()
      {
        return "<simuPOP::generalized step-wise mutator>" ;
      }

    private:
      /// probability to increase allele state
      double m_incProb;

      /// parameter for geometric gsm
      double m_p;

      /// the function to return random number
      PyObject* m_func;

  };

  /// mixed mutation model . has not been implemented.
  template<class Pop>
    class PyMutator: public Mutator<Pop>
  {
    public:
      PyMutator(vectorf rate=vectorf(), vectori atLoci=vectori(), UINT maxAllele=0,
        PyObject* func=NULL,
        string output=">", string outputExpr="", 
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t")
        :Mutator<Pop>( rate,atLoci, maxAllele,
        output, outputExpr, stage, begin, end, step, at, rep, grp, sep), m_func(NULL)
      {
        if( !PyCallable_Check(func))
          throw ValueError("Passed variable is not a callable python function.");

        Py_XINCREF(func);
        m_func = func;
      }

      ~PyMutator()
      {
        if( m_func != NULL )
          Py_DECREF(m_func);
      }

      /// CPPONLY
      PyMutator(const PyMutator& rhs):
      Mutator<Pop>(rhs), m_func(rhs.m_func)
      {
        if( m_func != NULL )
          Py_INCREF(m_func);
      }

      /// this function is very important
      virtual Operator<Pop>* clone() const
      {
        return new PyMutator<Pop>(*this);
      }

      virtual void mutate(Allele* allele)
      {
        PyObject* arglist = Py_BuildValue("(i)", *allele);
        PyObject* result = PyEval_CallObject(m_func, arglist);
        Py_DECREF(arglist);

        if( result == NULL)
        {
          PyErr_Print();
          throw ValueError("Function call failed.");
        }

        int resInt;
        PyObj_As_Int(result, resInt);
        DBG_DO(DBG_MUTATOR, cout << "Mutate " << static_cast<int>(*allele)
          << " to " << resInt << endl);
        *allele = static_cast<Allele>(resInt);
        Py_DECREF(result);
      }

      virtual string __repr__()
      {
        return "<simuPOP::python mutator>" ;
      }

    private:

      PyObject* m_func;
  };

  /// \brief point mutator
  /** mutate specified individuals at specified loci to spcified allele.
  I.e., this is a non-random mutator used to introduce disease etc.
  */
  template<class Pop>
    class PointMutator: public Operator<Pop>
  {
    public:
      /** \brief mutate once

      \param atLoci a vector of loci index.
      \param inds mutate 'inds' individuals
      \param toAllele mutate to 'toAllele'
      */
      PointMutator(
        vectori atLoci, Allele toAllele, vectorlu inds=vectorlu(),
        string output=">", string outputExpr="", 
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL, string sep="\t")
        :Operator<Pop>(output, outputExpr, stage, begin, end, step, at, rep, grp, sep),
        m_atLoci(atLoci), m_toAllele(toAllele), m_inds(inds)
      {

      }

      /// destructor
      virtual ~PointMutator()
      {
      }

      /// this function is very important
      virtual Operator<Pop>* clone() const
      {
        return new PointMutator<Pop>(*this);
      }

      /// apply!
      virtual bool apply(Pop& pop)
      {
        // mutate each mutable locus
        for( size_t i=0, iEnd=m_atLoci.size(); i < iEnd; ++i)
        {
          for( vectorlu::iterator ind = m_inds.begin();
            ind != m_inds.end(); ++ind)
          {
            *(pop.individual(*ind).genoBegin()+m_atLoci[i])= m_toAllele;
          }
        }                                         // each applicable loci

        return true;
      }

      virtual string __repr__()
      {
        return "<simuPOP::point mutator>" ;
      }

    private:

      /// applicable loci.
      vectori m_atLoci;
      Allele m_toAllele;
      vectorlu m_inds;
  };

}
#endif
