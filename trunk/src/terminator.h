/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$                                                      *
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

#ifndef _TERMINATOR_H
#define _TERMINATOR_H
/**
\file
\brief head file of class terminator: public Operator
*/
#include "population.h"
#include "operator.h"

#include <iostream>
#include <iomanip>

namespace simuPOP
{

  class terminator: public Operator
  {

    public:
      /// constructor. default to be always active.
      terminator(string message = "", string output=">", string outputExpr="",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL):
      Operator(output, outputExpr, stage, begin, end, step, at, rep, grp),
        m_message(message)
      {
      };

      /// destructor
      virtual ~terminator(){};

      virtual Operator* clone() const
      {
        return new terminator(*this);
      }

      string message()

      {
        return m_message;
      }

    private:
      /// message to print when terminated
      string m_message;
  };

  /// terminate according to a condition
  /// which can be, e.g.
  ///    any(alleleNum0) == 0
  ///    all(alleleNum1) > 0.5
  ///    alleleNum0{2} == 0
  /// etc.
  ///
  /// When the condition is true, a shared variable var="terminate" will be
  /// set to current generation.

  class terminateIf: public terminator
  {

    public:
      terminateIf(string condition="", string message="", string var="terminate",
        string output="", string outputExpr="",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL):
      terminator(message, output, outputExpr, stage, begin, end, step, at,
        rep, grp), m_expr(condition ), m_var(var)
      {
      }

      virtual Operator* clone() const
      {
        return new terminateIf(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::terminateIf>";
      }

      /// check all alleles in vector allele if they are fixed.
      virtual bool apply(population& pop)
      {
        // experssion return true
        m_expr.setLocalDict(pop.dict());

        if( m_expr.valueAsBool() == true )
        {
          // store the generations this replicate stops
          int gen = pop.gen();                    // mainVars().getVarAsInt("gen");

          pop.setIntVar(m_var, gen);
          if( !this->noOutput() )
          {
            ostream & out = this->getOstream(pop.dict());
            out << gen << endl;
            this->closeOstream();
          }
          if(this->message() != "")
            cout << this->message() << endl;
          return false;                           // return false, this replicate will be stopped
        }
        else
          return true;
      }

      virtual ~terminateIf(){};

    private:
      /// alleles to check. If empty, check all alleles.
      Expression m_expr;

      /// variable to set when terminated
      string m_var;
  };

  /// terminate according to a condition
  /// which can be, e.g.
  ///    any(alleleNum0) == 0
  ///    all(alleleNum1) > 0.5
  ///    alleleNum0{2} == 0
  /// etc.
  ///
  /// When the condition is true, a shared variable var="terminate" will be
  /// set to current generation.

  class continueIf: public terminator
  {

    public:
      continueIf(string condition="", string message="", string var="terminate",
        string output="", string outputExpr="",
        int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
        int rep=REP_ALL, int grp=GRP_ALL):
      terminator(message, output, outputExpr, stage, begin, end, step, at,
        rep, grp), m_expr(condition ), m_var(var)
      {
      }

      virtual Operator* clone() const
      {
        return new continueIf(*this);
      }

      virtual string __repr__()
      {
        return "<simuPOP::terminateIf>";
      }

      /// check all alleles in vector allele if they are fixed.
      virtual bool apply(population& pop)
      {
        // experssion return true
        m_expr.setLocalDict(pop.dict());

        if( m_expr.valueAsBool() == false )
        {
          // store the generations this replicate stops
          int gen = pop.gen();                    // mainVars().getVarAsInt("gen");

          pop.setIntVar(m_var, gen);
          if( !this->noOutput() )
          {
            ostream & out = this->getOstream(pop.dict());
            out << gen << endl;
            this->closeOstream();
          }
          if(this->message() != "")
            cout << this->message() << endl;
          return false;                           // return false, this replicate will be stopped
        }
        else
          return true;
      }

      virtual ~continueIf(){};

    private:
      /// alleles to check. If empty, check all alleles.
      Expression m_expr;

      /// variable to set when terminated
      string m_var;
  };

}
#endif
