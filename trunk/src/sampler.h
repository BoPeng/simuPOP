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

#ifndef _SAMPLER_H
#define _SAMPLER_H
/**
   \file
   \brief head file of class selector:public baseOperator
 */
#include "utility.h"
#include "operator.h"

#include "boost/tuple/tuple.hpp"
#include <numeric>
using std::min;

const string ASC_AS_Fields[2] = { "father_idx", "mother_idx" };

namespace simuPOP {
//
// /// draw a large pedigree sample
// class largePedigreeSample : public sample
// {
//
// public:
//  /// draw a large pedigree sample
//  /**
//   \param minTotalSize the minimum number of individuals in the sample
//   \param maxOffspring the maximum number of offspring a parent may have
//   \param minPedSize the minimal pedigree size. Default to \c 5.
//   \param minAffected the minimal number of affected individuals in each pedigree. Default to \c 0.
//   \param countOnly set variables about the number of affected sibpairs,
//      do not actually draw the sample.
//
//     Please refer to class \c sample for other parameter descriptions.
//   */
//  largePedigreeSample(vectoru size = vectoru(),
//                      unsigned minTotalSize = 0,
//                      unsigned maxOffspring = 5,
//                      unsigned minPedSize = 5,
//                      unsigned minAffected = 0,
//                      bool countOnly = false,
//                      const string & name = "sample", const string & nameExpr = "", UINT times = 1,
//                      const string & saveAs = "", const string & saveAsExpr = "",
//                      const string & format = "auto",
//                      int stage = PostMating, int begin = 0, int end = -1,
//                      int step = 1, vectorl at = vectorl(),
//                      const repList & rep = repList(), const subPopList & subPop = subPopList(),
//                      const vectorstr & infoFields = vectorstr (ASC_AS_Fields, ASC_AS_Fields + 2))
//      : sample(name, nameExpr, times, saveAs, saveAsExpr, format,
//               stage, begin, end, step, at, rep, subPop, infoFields),
//      m_size(size), m_minTotalSize(minTotalSize), m_maxOffspring(maxOffspring),
//      m_minPedSize(minPedSize), m_minAffected(minAffected),
//      m_countOnly(countOnly), m_validPedigrees()
//  {
//  }
//
//
//  /// destructor
//  virtual ~largePedigreeSample()
//  {
//  };
//
//  /// deep copy of a \c largePedigreeSample operator
//  virtual baseOperator * clone() const
//  {
//      return new largePedigreeSample(*this);
//  }
//
//
//  /// preparation before drawing a sample
//  virtual bool prepareSample(population & pop);
//
//  /// draw a a large pedigree sample
//  virtual population & drawsample(population & pop);
//
//  /// used by Python print function to print out the general information of the \c largePedigreeSample operator
//  virtual string __repr__()
//  {
//      return "<simuPOP::affected sibpair sample>" ;
//  }
//
//
// private:
//  typedef vector<boost::tuple<double, int> > pedArray;
//
//  /// sample size
//  vectoru m_size;
//
//  /// control total size
//  unsigned m_minTotalSize;
//
//  ///
//  unsigned m_maxOffspring;
//
//  ///
//  unsigned m_minPedSize;
//
//  ///
//  unsigned m_minAffected;
//
//  // do not draw sample
//  bool m_countOnly;
//
//  /// sibs for all subpopulations
//  /// we need to also save size information.
//  vector<pedArray> m_validPedigrees;
// };
//
// /// draw a nuclear family sample
// class nuclearFamilySample : public sample
// {
//
// public:
//  /// draw a nuclear family sample
//  /**
//     Please refer to class \c sample for parameter descriptions.
//   */
//  nuclearFamilySample(vectoru size = vectoru(),
//                      unsigned minTotalSize = 0,
//                      unsigned maxOffspring = 5,
//                      unsigned minPedSize = 5,
//                      unsigned minAffected = 0,
//                      bool countOnly = false,
//                      const string & name = "sample", const string & nameExpr = "", UINT times = 1,
//                      const string & saveAs = "", const string & saveAsExpr = "",
//                      const string & format = "auto",
//                      int stage = PostMating, int begin = 0, int end = -1,
//                      int step = 1, vectorl at = vectorl(),
//                      const repList & rep = repList(), const subPopList & subPop = subPopList(),
//                      const vectorstr & infoFields = vectorstr (ASC_AS_Fields, ASC_AS_Fields + 2))
//      : sample(name, nameExpr, times, saveAs, saveAsExpr, format,
//               stage, begin, end, step, at, rep, subPop, infoFields),
//      m_size(size), m_minTotalSize(minTotalSize), m_maxOffspring(maxOffspring),
//      m_minPedSize(minPedSize), m_minAffected(minAffected),
//      m_countOnly(countOnly), m_validPedigrees()
//  {
//  }
//
//
//  /// destructor
//  virtual ~nuclearFamilySample()
//  {
//  };
//
//  /// deep copy of a \c nuclearFamilySample operator
//  virtual baseOperator * clone() const
//  {
//      return new nuclearFamilySample(*this);
//  }
//
//
//  /// preparation before drawing a sample
//  virtual bool prepareSample(population & pop);
//
//  /// draw a nuclear family sample
//  virtual population & drawsample(population & pop);
//
//  /// used by Python print function to print out the general information of the \c nuclearFamilySample operator
//  virtual string __repr__()
//  {
//      return "<simuPOP::affected sibpair sample>" ;
//  }
//
//
// private:
//  typedef vector<boost::tuple<double, int> > pedArray;
//
//  /// sample size
//  vectoru m_size;
//
//  /// control total size
//  unsigned m_minTotalSize;
//
//  ///
//  unsigned m_maxOffspring;
//
//  ///
//  unsigned m_minPedSize;
//
//  ///
//  unsigned m_minAffected;
//
//  // do not draw sample
//  bool m_countOnly;
//
//  /// sibs for all subpopulations
//  /// we need to also save size information.
//  vector<pedArray> m_validPedigrees;
// };
//
// /// Python sampler
// /**
//    A Python sampler that generate a sample with given individuals.
//    This sampler accepts a Python array with elements that will be assigned to
//    individuals as their subpopulation IDs. Individuals with positive subpopulation IDs
//    will then be picked out and form a sample.
//
//    <funcForm>PySample</funcForm>
//  */
// class pySample : public sample
// {
//
// public:
//  /// create a Python sampler
//  /**
//   \param keep subpopulation IDs of all individuals
//   \param keepAncestralPop the number of ancestral populations that will be kept. If \c -1 is given,
//      keep all ancestral populations (default). If \c 0 is given, no ancestral population will be kept.
//
//     Please refer to class \c sample for other parameter descriptions.
//   */
//  pySample(PyObject * keep, int keepAncestralPops = -1,
//           const string & name = "sample", const string & nameExpr = "", UINT times = 1,
//           const string & saveAs = "", const string & saveAsExpr = "",   const string & format = "auto",
//           int stage = PostMating, int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
//           const repList & rep = repList(), const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr())
//      : sample(name, nameExpr, times, saveAs, saveAsExpr, format,
//               stage, begin, end, step, at, rep, subPop, infoFields),
//      m_keepAncestralPops(keepAncestralPops)
//  {
//      DBG_ASSERT(PyObj_Is_IntNumArray(keep), ValueError,
//          "Passed vector is not a Python/Numeric int array");
//      Py_INCREF(keep);
//      m_keep = keep;
//  }
//
//
//  /// destructor
//  virtual ~pySample()
//  {
//      if (m_keep != NULL)
//          Py_DECREF(m_keep);
//  }
//
//
//  /// CPPONLY
//  pySample(const pySample & rhs) :
//      sample(rhs),
//      m_keep(rhs.m_keep)
//  {
//      if (m_keep != NULL)
//          Py_INCREF(m_keep);
//  }
//
//
//  /// deep copy of a Python sampler
//  virtual baseOperator * clone() const
//  {
//      return new pySample(*this);
//  }
//
//
//  /// draw a Python sample
//  virtual population & drawsample(population & pop)
//  {
//      DBG_ASSERT(NumArray_Size(m_keep) >= static_cast<int>(pop.popSize()),
//          ValueError, "Given subpopid array has a length of "
//          + toStr(NumArray_Size(m_keep)) + " which is less than population size "
//          + toStr(pop.popSize()));
//
//      long * id = reinterpret_cast<long *>(NumArray_Data(m_keep));
//
//      for (size_t i = 0, iEnd = pop.popSize(); i < iEnd; ++i) {
//          DBG_ASSERT(static_cast<size_t>(id[i]) <= MaxSubPopID, ValueError,
//              "Subpop id exceeding maximum allowed subpopulations");
//          // convert from int to signed short
//          pop.ind(i).setSubPopID(static_cast<SubPopID>(id[i]));
//      }
//
//      return pop.newPopByIndID(m_keepAncestralPops);
//  }
//
//
//  /// used by Python print function to print out the general information of the Python sampler
//  virtual string __repr__()
//  {
//      return "<simuPOP::pySubset>" ;
//  }
//
//
// private:
//  PyObject * m_keep;
//
//  int m_keepAncestralPops;
// };
//
}
#endif
