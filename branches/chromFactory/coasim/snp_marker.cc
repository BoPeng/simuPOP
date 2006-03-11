/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "snp_marker.hh"

#ifndef CORE__DIST_FUNCTIONS_HH_INCLUDED
# include "dist_funcs.hh"
#endif
#ifndef CORE__RETIRED_INTERVAL_HH_INCLUDED
# include "retired_interval.hh"
#endif
#ifndef CORE__NODE_HH_INCLUDED
# include "node.hh"
#endif

namespace core
{

  Marker * SNPMarker::copy() const
  {
    return new SNPMarker(*this);
  }

  bool SNPMarker::run_first() const
  {
    return false;
  }

  int SNPMarker::default_value() const
  {
    return 0;
  }

  namespace
  {
    class SNPMutator : public Mutator
    {

      const SNPMarker &m_marker;

      // allowed frequencies, translated into leaf node counts
      unsigned int m_low_leaf_count;
      unsigned int m_high_leaf_count;

      double m_mutation_point;                    // the point on the "surface" where
      // the mutation sits
      double m_surface_so_far;                    // the surface seen so far

      public:
        SNPMutator(const SNPMarker &marker,
          unsigned int low_leaf_count,
          unsigned int high_leaf_count,
          double mutation_point);

        bool edge_has_mutation(double parent_time, double child_time);
        int  mutate_to(const Node &n, int parent_allele);
        virtual int  mutate(const Node &parent,
          const Node &child,
          int parent_allele);
    };

    SNPMutator::SNPMutator(const SNPMarker &marker,
      unsigned int low_leaf_count,
      unsigned int high_leaf_count,
      double mutation_point)
      : m_marker(marker),
      m_low_leaf_count(low_leaf_count), m_high_leaf_count(high_leaf_count),
      m_mutation_point(mutation_point), m_surface_so_far(0.0)
    {
    }

    bool SNPMutator::edge_has_mutation(double parent_time, double child_time)
    {
      double edge_length = parent_time - child_time;
      bool mutate = (m_surface_so_far <= m_mutation_point
        and
        m_mutation_point < m_surface_so_far+edge_length);
      m_surface_so_far += edge_length;
      return mutate;
    }

    int SNPMutator::mutate_to(const Node &child, int parent_allele)
    {
      // check frequency
      unsigned int leaf_count = child.leaves_at_point(m_marker.position());
      if ((leaf_count < m_low_leaf_count) or (m_high_leaf_count < leaf_count))
        throw Mutator::retry_mutation();
      return !parent_allele;
    }

    int SNPMutator::mutate(const Node &parent, const Node &child,
      int parent_allele)
    {
      if (edge_has_mutation(parent.time(), child.time()))
        return mutate_to(child, parent_allele);
      else
        return parent_allele;
    }
  }

  static inline void swap(unsigned int &i, unsigned int &j)
    { unsigned int tmp = i; i = j; j = tmp; }

  Mutator *
    SNPMarker::create_mutator(const Configuration   &conf,
    const RetiredInterval &ri) const
  {
    unsigned int low_leaf_count
      = static_cast<unsigned int>(ceil(m_low_freq*conf.no_leaves()));
    unsigned int high_leaf_count
      = static_cast<unsigned int>(floor(m_high_freq*conf.no_leaves()));

    // possible due to ceil/floor
    if (high_leaf_count < low_leaf_count) swap(low_leaf_count,high_leaf_count);

    double mutation_point = ri.surface() * Distribution_functions::uniform();

    return new SNPMutator(*this, low_leaf_count, high_leaf_count, mutation_point);
  }

  const char * SNPMarker::type() const
  {
    return "snp";
  }
}
